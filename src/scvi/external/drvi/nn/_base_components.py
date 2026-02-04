from __future__ import annotations

import collections
import math
from typing import TYPE_CHECKING

import torch
from torch import nn
from torch.distributions import Normal
from torch.nn import functional as F

from scvi.external.drvi.nn_modules.embedding import MultiEmbedding
from scvi.external.drvi.nn_modules.freezable import FreezableBatchNorm1d, FreezableLayerNorm
from scvi.external.drvi.nn_modules.gradients import GradientScaler
from scvi.external.drvi.nn_modules.layer.factory import FCLayerFactory, LayerFactory
from scvi.external.drvi.nn_modules.layer.linear_layer import StackedLinearLayer
from scvi.external.drvi.nn_modules.noise_model import NoiseModel

if TYPE_CHECKING:
    from collections.abc import Callable, Iterable, Sequence
    from typing import Any, Literal


def _identity(x: torch.Tensor) -> torch.Tensor:
    return x


class FCLayers(nn.Module):
    """A helper class to build fully-connected layers for a neural network.

    Parameters
    ----------
    layers_dim
        Number of nodes in layers including input and output dimensions.
    n_cat_list
        A list containing, for each category of interest,
        the number of categories. Each category will be
        included using a one-hot encoding.
    dropout_rate
        Dropout rate to apply to each of the hidden layers.
    n_split
        The size of split if input is a 3d tensor otherwise 1.
        This parameter is required to handle batch normalization.
    reuse_weights
        Whether to reuse weights when having multiple splits (defined per layer). None if no reuse / split.
    use_batch_norm
        Whether to have `BatchNorm` layers or not.
    affine_batch_norm
        Whether to have affine transformation in `BatchNorm` layers.
    use_layer_norm
        Whether to have `LayerNorm` layers or not.
    use_activation
        Whether to have layer activation or not.
    bias
        Whether to learn bias in linear layers or not.
    inject_covariates
        Whether to inject covariates in each layer, or just the first.
    activation_fn
        Which activation function to use.
    layer_factory
        A layer Factory instance to build projection layers based on.
    layers_location
        An indicator to tell the class where in the architecture these layers reside.
    covariate_modeling_strategy
        The strategy model to consider covariates.
    covariate_embs_dim
        Dimensions for covariate embeddings when using embedding strategies.
    """

    def __init__(
        self,
        layers_dim: Sequence[int],
        n_cat_list: Iterable[int] | None = None,
        dropout_rate: float = 0.1,
        n_split: int = 1,
        reuse_weights: Sequence[bool] | None = None,
        use_batch_norm: bool = True,
        affine_batch_norm: bool = True,
        use_layer_norm: bool = False,
        use_activation: bool = True,
        bias: bool = True,
        inject_covariates: bool = True,
        activation_fn: type[nn.Module] = nn.ELU,
        layer_factory: LayerFactory | None = None,
        layers_location: Literal["intermediate", "first", "last"] = "intermediate",
        covariate_modeling_strategy: Literal[
            "one_hot",
            "emb",
            "emb_shared",
            "one_hot_linear",
            "emb_linear",
            "emb_shared_linear",
        ] = "one_hot",
        covariate_embs_dim: Iterable[int] = (),
    ) -> None:
        super().__init__()
        assert n_split > 0

        self.inject_covariates = inject_covariates
        if covariate_modeling_strategy.endswith("_linear"):
            self.covariate_projection_modeling = "linear"
            self.covariate_vector_modeling = covariate_modeling_strategy[: -len("_linear")]
        else:
            self.covariate_projection_modeling = "cat"
            self.covariate_vector_modeling = covariate_modeling_strategy
        layer_factory = layer_factory or FCLayerFactory()

        self.n_cat_list = list(n_cat_list) if n_cat_list is not None else []
        if self.covariate_vector_modeling == "one_hot":
            covariate_embs_dim = [n_cat if n_cat > 1 else 0 for n_cat in self.n_cat_list]
        else:
            covariate_embs_dim = list(covariate_embs_dim)
            assert len(covariate_embs_dim) == len(self.n_cat_list)

        self.injectable_layers = []
        self.linear_projections = []
        self.stacked_layers = []
        self.linear_batch_projections = []

        def is_intermediate(i: int) -> bool:
            assert layers_location in ["intermediate", "first", "last"]
            if layers_location == "first" and i == 0:
                return False
            if layers_location == "last" and i == len(layers_dim) - 2:
                return False
            return True

        def inject_into_layer(layer_num: int) -> bool:
            user_cond = layer_num == 0 or (layer_num > 0 and self.inject_covariates)
            return user_cond

        def get_projection_layer(n_in: int, n_out: int, i: int) -> list[nn.Module]:
            output = []
            layer_needs_injection = False
            if (reuse_weights is not None) and (not all(reuse_weights)):
                assert n_split > 1
                if self.covariate_projection_modeling not in ["cat", "linear"]:
                    raise NotImplementedError()

            if len(self.n_cat_list) > 0 and inject_into_layer(i):
                layer_needs_injection = True
                cat_dim = sum(covariate_embs_dim)
                if self.covariate_vector_modeling == "emb":
                    batch_emb = MultiEmbedding(
                        self.n_cat_list, covariate_embs_dim, init_method="normal", max_norm=1.0
                    )
                    output.append(batch_emb)
                if self.covariate_projection_modeling == "cat":
                    n_in += cat_dim
                elif self.covariate_projection_modeling == "linear":
                    if (reuse_weights is None) or reuse_weights[i]:
                        linear_batch_projection = nn.Linear(cat_dim, n_out, bias=False)
                    else:
                        linear_batch_projection = StackedLinearLayer(
                            n_split, cat_dim, n_out, bias=False
                        )
                    output.append(linear_batch_projection)
                    self.linear_batch_projections.append(linear_batch_projection)
                else:
                    raise NotImplementedError()
            if (reuse_weights is None) or reuse_weights[i]:
                layer = layer_factory.get_normal_layer(
                    n_in,
                    n_out,
                    bias=bias,
                    intermediate_layer=is_intermediate(i),
                )
            else:
                layer = layer_factory.get_stacked_layer(
                    n_split,
                    n_in,
                    n_out,
                    bias=bias,
                    intermediate_layer=is_intermediate(i),
                )
                self.stacked_layers.append(layer)
            self.linear_projections.append(layer)
            if layer_needs_injection:
                self.injectable_layers.append(layer)
            output.append(layer)
            return output

        def get_normalization_layers(n_out: int) -> list[nn.Module]:
            output = []
            if n_split == 1:
                if use_batch_norm:
                    # non-default params come from defaults in original Tensorflow implementation
                    output.append(
                        FreezableBatchNorm1d(
                            n_out, momentum=0.01, eps=0.001, affine=affine_batch_norm
                        )
                    )
                if use_layer_norm:
                    output.append(FreezableLayerNorm(n_out, elementwise_affine=False))
            else:
                if use_batch_norm:
                    # non-default params come from defaults in original Tensorflow implementation
                    output.append(
                        FreezableBatchNorm1d(
                            n_out * n_split, momentum=0.01, eps=0.001, affine=affine_batch_norm
                        )
                    )
                if use_layer_norm:
                    output.append(FreezableLayerNorm(n_out, elementwise_affine=False))
                    # The following logic is wrong
                    # output.append(FreezableLayerNorm([n_split, n_out], elementwise_affine=False))
            return output

        self.fc_layers = nn.Sequential(
            collections.OrderedDict(
                [
                    (
                        f"Layer {i}",
                        nn.ModuleList(
                            [
                                p
                                for p in [
                                    *get_projection_layer(n_in, n_out, i),
                                    *get_normalization_layers(n_out),
                                    activation_fn() if use_activation else None,
                                    nn.Dropout(p=dropout_rate) if dropout_rate > 0 else None,
                                ]
                                if p is not None
                            ]
                        ),
                    )
                    for i, (n_in, n_out) in enumerate(
                        zip(layers_dim[:-1], layers_dim[1:], strict=True)
                    )
                ]
            )
        )

    def set_online_update_hooks(
        self, previous_n_cats_per_cov: Sequence[int], n_cats_per_cov: Sequence[int]
    ) -> None:
        """Set online update hooks for handling new categories.

        Parameters
        ----------
        previous_n_cats_per_cov
            Previous number of categories per covariate.
        n_cats_per_cov
            New number of categories per covariate.
        """
        if sum(previous_n_cats_per_cov) == sum(n_cats_per_cov):
            print("Nothing to make hook for!")
            return

        def make_hook_function(weight: torch.Tensor) -> Callable[[torch.Tensor], torch.Tensor]:
            w_size = weight.size()
            sum_n_cats_per_cov = sum([n_cat if n_cat > 1 else 0 for n_cat in n_cats_per_cov])
            if weight.dim() == 2:
                # 2D tensors
                with torch.no_grad():
                    # Freeze gradients for normal nodes
                    if w_size[1] == sum_n_cats_per_cov:
                        transfer_mask = []
                    else:
                        transfer_mask = [
                            torch.zeros(
                                [w_size[0], w_size[1] - sum_n_cats_per_cov], device=weight.device
                            )
                        ]
                    # Iterate over the categories and Freeze old caterogies and make new ones trainable
                    for n_cat_new, n_cat_old in zip(
                        n_cats_per_cov, previous_n_cats_per_cov, strict=False
                    ):
                        n_cat_new = n_cat_new if n_cat_new > 1 else 0
                        n_cat_old = n_cat_old if n_cat_old > 1 else 0
                        transfer_mask.append(
                            torch.zeros([w_size[0], n_cat_old], device=weight.device)
                        )
                        if n_cat_new > n_cat_old:
                            transfer_mask.append(
                                torch.ones(
                                    [w_size[0], n_cat_new - n_cat_old], device=weight.device
                                )
                            )
                    transfer_mask = torch.cat(transfer_mask, dim=1)
            elif weight.dim() == 3:
                # 3D tensors
                with torch.no_grad():
                    # Freeze gradients for normal nodes
                    if w_size[1] == sum_n_cats_per_cov:
                        transfer_mask = []
                    else:
                        transfer_mask = [
                            torch.zeros(
                                [w_size[0], w_size[1] - sum_n_cats_per_cov, w_size[2]],
                                device=weight.device,
                            )
                        ]
                    # Iterate over the categories and Freeze old caterogies and make new ones trainable
                    for n_cat_new, n_cat_old in zip(
                        n_cats_per_cov, previous_n_cats_per_cov, strict=False
                    ):
                        n_cat_new = n_cat_new if n_cat_new > 1 else 0
                        n_cat_old = n_cat_old if n_cat_old > 1 else 0
                        transfer_mask.append(
                            torch.zeros([w_size[0], n_cat_old, w_size[2]], device=weight.device)
                        )
                        if n_cat_new > n_cat_old:
                            transfer_mask.append(
                                torch.ones(
                                    [w_size[0], n_cat_new - n_cat_old, w_size[2]],
                                    device=weight.device,
                                )
                            )
                    transfer_mask = torch.cat(transfer_mask, dim=1)
            else:
                raise NotImplementedError()

            def _hook_fn_injectable(grad: torch.Tensor) -> torch.Tensor:
                return grad * transfer_mask

            return _hook_fn_injectable

        for layers in self.fc_layers:
            for layer in layers:
                if self.covariate_vector_modeling == "emb_shared":
                    # Nothing to do here :)
                    pass
                elif self.covariate_vector_modeling == "emb":
                    # Freeze everything but embs (new embs)
                    if isinstance(layer, MultiEmbedding):
                        assert tuple(layer.num_embeddings) == tuple(n_cats_per_cov)
                        print(f"Freezing old categories for {layer}")
                        layer.freeze_top_embs(previous_n_cats_per_cov)
                    # No need to handle others. For them required_grad is already set to False
                elif self.covariate_vector_modeling == "one_hot":
                    assert self.covariate_projection_modeling in ["cat", "linear"]
                    # Freeze everything but linears right after one_hot (new weights)
                    if (
                        self.covariate_projection_modeling == "cat"
                        and layer in self.injectable_layers
                    ) or (
                        self.covariate_projection_modeling == "linear"
                        and layer in self.linear_batch_projections
                    ):
                        assert layer.weight.requires_grad
                        print(
                            f"Registering backward hook parameter with shape {layer.weight.size()}"
                        )
                        layer.weight.register_hook(make_hook_function(layer.weight))
                        if layer.bias is not None:
                            assert not layer.bias.requires_grad
                else:
                    raise NotImplementedError()

    def forward(
        self,
        x: torch.Tensor,
        cat_full_tensor: torch.Tensor | None,
        output_subset_indices: torch.Tensor | None = None,
    ) -> torch.Tensor:
        """Forward computation on ``x``.

        Parameters
        ----------
        x
            Tensor of values with shape ``(batch_size, n_in,)`` or ``(batch_size, n_split, n_in)``.
        cat_full_tensor
            Tensor of category membership(s) for this sample.
        output_subset_indices
            The subset of output dims.

        Returns
        -------
        torch.Tensor
            Tensor of shape ``(batch_size, n_out,)`` or ``(batch_size, n_split, n_out)``.
        """
        if self.covariate_vector_modeling == "one_hot":
            concat_list = []
            if cat_full_tensor is not None:
                cat_list = torch.split(cat_full_tensor, 1, dim=1)
            else:
                cat_list = ()

            if len(self.n_cat_list) > len(cat_list):
                raise ValueError("nb. categorical args provided doesn't match init. params.")
            for n_cat, cat in zip(self.n_cat_list, cat_list, strict=False):
                if n_cat and cat is None:
                    raise ValueError("cat not provided while n_cat != 0 in init. params.")
                if n_cat > 1:  # n_cat = 1 will be ignored - no additional information
                    concat_list += [F.one_hot(cat.long().squeeze(-1), n_cat).float()]
        elif self.covariate_vector_modeling == "emb_shared":
            concat_list = [cat_full_tensor]
        else:
            concat_list = []

        def dimension_transformation(t: torch.Tensor) -> torch.Tensor:
            if x.dim() == t.dim():
                return t
            if x.dim() == 3 and t.dim() == 2:
                return t.unsqueeze(dim=1).expand(-1, x.shape[1], -1)
            raise NotImplementedError()

        for layers in self.fc_layers:
            concat_list_layer = concat_list
            projected_batch_layer = None
            for layer in layers:
                if layer is not None:
                    if isinstance(layer, nn.BatchNorm1d):
                        if x.dim() == 3:
                            x = layer(x.reshape(x.shape[0], -1)).reshape(x.shape)
                        else:
                            x = layer(x)
                    elif isinstance(layer, MultiEmbedding):
                        assert self.covariate_vector_modeling in ["emb"]
                        assert len(concat_list) == 0
                        concat_list_layer = [layer(cat_full_tensor.int())]
                    elif layer in self.linear_batch_projections:
                        assert self.covariate_projection_modeling == "linear"
                        projected_batch_layer = layer(torch.cat(concat_list_layer, dim=-1))
                    else:
                        kwargs = {}
                        if layer == self.linear_projections[-1]:
                            kwargs["output_subset"] = output_subset_indices
                        to_be_added = None
                        if layer in self.injectable_layers:
                            if self.covariate_projection_modeling == "cat":
                                if len(concat_list_layer) > 0:
                                    current_cat_tensor = dimension_transformation(
                                        torch.cat(concat_list_layer, dim=-1)
                                    )
                                    x = torch.cat((x, current_cat_tensor), dim=-1)
                            elif self.covariate_projection_modeling in ["linear"]:
                                to_be_added = dimension_transformation(projected_batch_layer)
                            else:
                                raise NotImplementedError()
                        x = layer(x, **kwargs)
                        if to_be_added is not None:
                            x = x + to_be_added
        return x


# Encoder
class Encoder(nn.Module):
    """Encode data of ``n_input`` dimensions into a latent space of ``n_output`` dimensions.

    Uses a fully-connected neural network of ``n_hidden`` layers.

    Parameters
    ----------
    n_input
        The dimensionality of the input (data space).
    n_output
        The dimensionality of the output (latent space).
    layers_dim
        The number of nodes per hidden layer as a sequence.
    n_cat_list
        A list containing the number of categories
        for each category of interest. Each category will be
        included using a one-hot encoding.
    n_continuous_cov
        The number of continuous covariates.
    inject_covariates
        Whether to inject covariates in each layer, or just the first.
    use_batch_norm
        Whether to use batch norm in layers.
    affine_batch_norm
        Whether to use affine in batch norms.
    use_layer_norm
        Whether to use layer norm in layers.
    input_dropout_rate
        Dropout rate to apply to the input.
    dropout_rate
        Dropout rate to apply to each of the hidden layers.
    distribution
        Distribution of z.
    var_eps
        Minimum value for the variance; used for numerical stability.
    var_activation
        Callable used to ensure positivity of the variance.
    mean_activation
        Callable used to apply activation to the mean.
    layer_factory
        A layer Factory instance to build projection layers based on.
    covariate_modeling_strategy
        The strategy model to consider covariates.
    categorical_covariate_dims
        Dimensions for covariate embeddings when using embedding strategies.
    return_dist
        Whether to return the distribution or just the parameters.
    **kwargs
        Additional keyword arguments.
    """

    def __init__(
        self,
        n_input: int,
        n_output: int,
        layers_dim: Sequence[int] = (128,),
        n_cat_list: Iterable[int] | None = None,
        n_continuous_cov: int = 0,
        inject_covariates: bool = True,
        use_batch_norm: bool = True,
        affine_batch_norm: bool = True,
        use_layer_norm: bool = False,
        input_dropout_rate: float = 0.0,
        dropout_rate: float = 0.1,
        distribution: str = "normal",
        var_eps: float = 1e-4,
        var_activation: Callable | Literal["exp", "pow2", "2sig"] = "exp",
        mean_activation: Callable | str = "identity",
        layer_factory: LayerFactory | None = None,
        covariate_modeling_strategy: Literal[
            "one_hot",
            "emb",
            "emb_shared",
            "one_hot_linear",
            "emb_linear",
            "emb_shared_linear",
        ] = "one_hot",
        categorical_covariate_dims: Sequence[int] = (),
        return_dist: bool = False,
        **kwargs,
    ) -> None:
        super().__init__()
        self.distribution = distribution
        self.var_eps = var_eps
        self.input_dropout = nn.Dropout(p=input_dropout_rate)
        self.fully_deterministic = False

        all_layers_dim = [n_input + n_continuous_cov] + list(layers_dim) + [n_output]
        if len(layers_dim) >= 1:
            self.encoder = FCLayers(
                layers_dim=all_layers_dim[:-1],
                n_cat_list=n_cat_list,
                dropout_rate=dropout_rate,
                use_batch_norm=use_batch_norm,
                affine_batch_norm=affine_batch_norm,
                use_layer_norm=use_layer_norm,
                inject_covariates=inject_covariates,
                layer_factory=layer_factory,
                layers_location="first",
                covariate_modeling_strategy=covariate_modeling_strategy,
                covariate_embs_dim=categorical_covariate_dims,
                **kwargs,
            )
        else:
            self.register_parameter("encoder", None)
            inject_covariates = True
        self.mean_encoder = FCLayers(
            layers_dim=all_layers_dim[-2:],
            n_cat_list=n_cat_list if inject_covariates else [],
            use_activation=False,
            use_batch_norm=False,
            use_layer_norm=False,
            bias=True,
            dropout_rate=0,
            layer_factory=layer_factory,
            layers_location="intermediate" if len(layers_dim) >= 1 else "first",
            covariate_modeling_strategy=covariate_modeling_strategy,
            covariate_embs_dim=categorical_covariate_dims if inject_covariates else [],
            **kwargs,
        )
        self.var_encoder = FCLayers(
            layers_dim=all_layers_dim[-2:],
            n_cat_list=n_cat_list if inject_covariates else [],
            use_activation=False,
            use_batch_norm=False,
            use_layer_norm=False,
            bias=True,
            dropout_rate=0,
            layer_factory=layer_factory,
            layers_location="intermediate" if len(layers_dim) >= 1 else "first",
            covariate_modeling_strategy=covariate_modeling_strategy,
            covariate_embs_dim=categorical_covariate_dims if inject_covariates else [],
            **kwargs,
        )
        self.return_dist = return_dist

        if distribution == "ln":
            self.z_transformation = nn.Softmax(dim=-1)
        else:
            self.z_transformation = _identity
        if var_activation == "exp":
            self.var_activation = torch.exp
        elif var_activation == "pow2":
            self.var_activation = lambda x: torch.pow(x, 2)
        elif var_activation == "2sig":
            self.var_activation = lambda x: 2 * torch.sigmoid(x)
        else:
            assert callable(var_activation)
            self.var_activation = var_activation

        if mean_activation == "identity":
            self.mean_activation = nn.Identity()
        elif mean_activation == "relu":
            self.mean_activation = nn.ReLU()
        elif mean_activation.startswith("leaky_relu"):
            if mean_activation == "leaky_relu":
                mean_activation = "leaky_relu_0.01"
            slope = float(mean_activation.split("leaky_relu_")[1])
            self.mean_activation = nn.LeakyReLU(negative_slope=slope)
        elif mean_activation.startswith("elu"):
            if mean_activation == "elu":
                mean_activation = "elu_1.0"
            alpha = float(mean_activation.split("elu_")[1])
            self.mean_activation = nn.ELU(alpha=alpha)
        elif mean_activation.startswith("celu"):
            if mean_activation == "celu":
                mean_activation = "celu_1.0"
            alpha = float(mean_activation.split("celu_")[1])
            self.mean_activation = nn.CELU(alpha=alpha)
        else:
            assert callable(mean_activation)
            self.mean_activation = mean_activation

    def _generate_variational_distribution(self, qz_m: torch.Tensor, qz_v: torch.Tensor) -> Normal:
        """Generate the variational distribution.

        Parameters
        ----------
        qz_m
            Mean of the variational distribution with shape ``(batch_size, n_latent)``.
        qz_v
            Variance of the variational distribution with shape ``(batch_size, n_latent)``.

        Returns
        -------
        Normal
            Variational distribution.
        """
        return Normal(qz_m, qz_v.sqrt())

    def _sample_from_variational(self, dist: Normal) -> torch.Tensor:
        """Sample from the variational distribution.

        Parameters
        ----------
        dist
            Variational distribution.

        Returns
        -------
        torch.Tensor
            Sampled latent tensor with shape ``(batch_size, n_latent)``.
        """
        if self.fully_deterministic:
            z = dist.mean
        else:
            z = dist.rsample()
        return self.z_transformation(z)

    def forward(
        self,
        x: torch.Tensor,
        cat_full_tensor: torch.Tensor | None,
        cont_full_tensor: torch.Tensor | None = None,
    ) -> tuple[torch.Tensor, torch.Tensor, torch.Tensor] | tuple[Normal, torch.Tensor]:
        r"""Forward computation on ``x``.

         #. Encodes the data into latent space using the encoder network
         #. Generates a mean \\( q_m \\) and variance \\( q_v \\)
         #. Samples a new value from an i.i.d. multivariate normal \\( \\sim Ne(q_m, \\mathbf{I}q_v) \\)

        Parameters
        ----------
        x
            Tensor with shape (batch_size, n_input).
        cat_full_tensor
            Tensor containing encoding of categorical variables of size n_batch x n_total_cat.
        cont_full_tensor
            Tensor containing continuous covariates.

        Returns
        -------
        tuple
            If return_dist=False: (mean, variance, sample) tensors of shape (batch_size, n_latent).
            If return_dist=True: (distribution, sample) where distribution is a Normal distribution.

        """
        if cont_full_tensor is not None:
            x = torch.cat((x, cont_full_tensor), dim=-1)

        # Parameters for latent distribution
        q = self.encoder(self.input_dropout(x), cat_full_tensor) if self.encoder is not None else x
        q_m = self.mean_activation(self.mean_encoder(q, cat_full_tensor))
        q_v = self.var_activation(self.var_encoder(q, cat_full_tensor)) + self.var_eps

        dist = self._generate_variational_distribution(q_m, q_v)
        if self.return_dist:
            return dist

        latent = self._sample_from_variational(dist)
        return q_m, q_v, latent


# Decoder
class DecoderDRVI(nn.Module):
    """Decodes data from latent space of ``n_input`` dimensions into ``n_output`` dimensions.

    Uses a fully-connected neural network of ``n_hidden`` layers.

    Parameters
    ----------
    n_input
        The dimensionality of the input (latent space).
    n_output
        The dimensionality of the output (data space).
    gene_likelihood_module
        Module defining the noise model for gene expression.
    n_cat_list
        A list containing the number of categories
        for each category of interest. Each category will be
        included using a one-hot encoding.
    n_continuous_cov
        The number of continuous covariates.
    n_split
        Number of splits in the latent space.
    split_aggregation
        How to aggregate splits in the last layer of the decoder.
    split_method
        How to make splits:
        - "split" : Split the latent space
        - "power" : Transform the latent space to n_split vectors of size n_latent
        - "split_map" : Split the latent space then map each to latent space using unique transformations
    reuse_weights
        Where to reuse the weights of the decoder layers when using splitting.
    layers_dim
        The number of nodes per hidden layer as a sequence.
    dropout_rate
        Dropout rate to apply to each of the hidden layers.
    inject_covariates
        Whether to inject covariates in each layer, or just the first.
    use_batch_norm
        Whether to use batch norm in layers.
    affine_batch_norm
        Whether to use affine in batch norms.
    use_layer_norm
        Whether to use layer norm in layers.
    layer_factory
        A layer Factory instance for building layers.
    covariate_modeling_strategy
        The strategy model takes to model covariates.
    categorical_covariate_dims
        Dimensions for categorical covariate embeddings.
    last_layer_gradient_scale
        Gradient scale for the last layer of the decoder.
    **kwargs
        Additional keyword arguments.
    """

    def __init__(
        self,
        n_input: int,
        n_output: int,
        gene_likelihood_module: NoiseModel,
        n_cat_list: Iterable[int] | None = None,
        n_continuous_cov: int = 0,
        n_split: int = 1,
        split_aggregation: Literal["sum", "logsumexp", "max"] = "logsumexp",
        split_method: Literal["split", "power", "split_map", "split_diag"] | str = "split",
        reuse_weights: Literal[
            "everywhere", "last", "intermediate", "nowhere", "not_first"
        ] = "everywhere",
        layers_dim: Sequence[int] = (128,),
        dropout_rate: float = 0.1,
        inject_covariates: bool = True,
        use_batch_norm: bool = False,
        affine_batch_norm: bool = True,
        use_layer_norm: bool = False,
        layer_factory: LayerFactory | None = None,
        covariate_modeling_strategy: Literal[
            "one_hot",
            "emb",
            "emb_shared",
            "one_hot_linear",
            "emb_linear",
            "emb_shared_linear",
        ] = "one_hot",
        categorical_covariate_dims: Sequence[int] = (),
        last_layer_gradient_scale: float = 1.0,
        **kwargs,
    ) -> None:
        super().__init__()
        self.n_output = n_output
        self.gene_likelihood_module = gene_likelihood_module
        self.n_split = n_split
        self.split_aggregation = split_aggregation
        self.split_method = split_method

        if self.n_split == -1 or self.n_split == 1:
            assert reuse_weights
            self.n_split = n_split = 1

        self.effective_dim = self._initialize_splitting(n_input)

        assert reuse_weights in ["everywhere", "last", "intermediate", "nowhere", "not_first"]
        if reuse_weights in ["everywhere", "intermediate"]:
            intermediate_layers_reuse_weights = [True] * len(layers_dim)
        elif reuse_weights == "not_first":
            intermediate_layers_reuse_weights = [True] * len(layers_dim)
            intermediate_layers_reuse_weights[0] = False
        else:  # last, nowhere
            intermediate_layers_reuse_weights = None
        last_layers_reuse_weights = reuse_weights in ["everywhere", "last", "not_first"]

        all_layers_dim = [self.effective_dim + n_continuous_cov] + list(layers_dim) + [n_output]
        if len(layers_dim) >= 1:
            self.px_shared_decoder = FCLayers(
                layers_dim=all_layers_dim[:-1],
                n_split=self.n_split,
                reuse_weights=intermediate_layers_reuse_weights,
                n_cat_list=n_cat_list,
                dropout_rate=dropout_rate,
                inject_covariates=inject_covariates,
                use_batch_norm=use_batch_norm,
                affine_batch_norm=affine_batch_norm,
                use_layer_norm=use_layer_norm,
                layer_factory=layer_factory,
                layers_location="intermediate",
                covariate_modeling_strategy=covariate_modeling_strategy,
                covariate_embs_dim=categorical_covariate_dims,
                **kwargs,
            )
        else:
            assert reuse_weights == "nowhere"
            self.register_parameter("px_shared_decoder", None)
            inject_covariates = True

        if last_layer_gradient_scale != 1.0:
            self.last_layer_gradient_scaler = GradientScaler(last_layer_gradient_scale)
        else:
            self.register_parameter("last_layer_gradient_scaler", None)

        params_for_likelihood = self.gene_likelihood_module.parameters
        params_nets = {}
        for param_name, param_info in params_for_likelihood.items():
            if param_info.startswith("fixed="):
                params_nets[param_name] = torch.nn.Parameter(
                    torch.tensor(float(param_info.split("=")[1])), requires_grad=False
                )
            elif param_info == "no_transformation":
                params_nets[param_name] = FCLayers(
                    layers_dim=all_layers_dim[-2:],
                    n_split=self.n_split,
                    reuse_weights=[last_layers_reuse_weights],
                    n_cat_list=n_cat_list if inject_covariates else [],
                    use_activation=False,
                    use_batch_norm=False,
                    use_layer_norm=False,
                    bias=True,
                    dropout_rate=0,
                    layer_factory=layer_factory,
                    layers_location="last",
                    covariate_modeling_strategy=covariate_modeling_strategy,
                    covariate_embs_dim=categorical_covariate_dims if inject_covariates else [],
                    **kwargs,
                )
            elif param_info == "per_feature":
                params_nets[param_name] = torch.nn.Parameter(torch.randn(n_output))
            else:
                raise NotImplementedError()
        self.params_nets = nn.ParameterDict(params_nets)

    def _initialize_splitting(self, n_input: int) -> None:
        """Initialize the splitting."""
        effective_dim = n_input
        if self.n_split == 1:
            pass
        elif self.split_method == "split":
            assert n_input % self.n_split == 0
            effective_dim = n_input // self.n_split
        elif self.split_method.startswith("split_map"):
            assert n_input % self.n_split == 0
            if "@" in self.split_method:
                effective_dim = int(self.split_method.split("@")[-1])
            effective_split_size = n_input // self.n_split
            self.split_transformation_weight = nn.Parameter(
                torch.randn(self.n_split, effective_split_size, effective_dim)
            )
            self.split_transformation_weight.data /= float(effective_split_size) ** 0.5
        elif self.split_method == "split_diag":
            assert n_input == self.n_split
        elif self.split_method.startswith("power"):
            if "@" in self.split_method:
                effective_dim = int(self.split_method.split("@")[-1])
            self.split_transformation = nn.Sequential(
                nn.Linear(n_input, effective_dim * self.n_split), nn.ReLU()
            )
        else:
            raise NotImplementedError()
        return effective_dim

    def _apply_split_transformation(self, z: torch.Tensor) -> torch.Tensor:
        """Apply split transformation to the latent tensor according to split method.

        Parameters
        ----------
        z
            Sampled latent tensor with shape ``(batch_size, n_latent)``.

        Returns
        -------
        torch.Tensor
            Transformed latent tensor with shape ``(batch_size, n_split, effective_dim)`` or ``(batch_size, effective_dim)``
        """
        batch_size = z.shape[0]

        if self.n_split > 1:
            if self.split_method == "split":
                z = torch.reshape(z, (batch_size, self.n_split, -1))
            elif self.split_method == "split_diag":
                z = torch.diag_embed(
                    z
                )  # (batch_size, n_latent) -> (batch_size, n_latent, n_latent)
            elif self.split_method.startswith("power"):
                z = self.split_transformation(z)
                z = torch.reshape(z, (batch_size, self.n_split, -1))
            elif self.split_method.startswith("split_map"):
                z = torch.reshape(z, (batch_size, self.n_split, -1))
                z = torch.einsum("bsd,sdn->bsn", z, self.split_transformation_weight)

        return z

    def _aggregate_split(self, x: torch.Tensor, dim: int = -2) -> torch.Tensor:
        """Aggregate the split.

        Parameters
        ----------
        x
            Tensor with shape ``(batch_size, n_split, n_features)``.
        dim
            Dimension to aggregate over.

        Returns
        -------
        torch.Tensor
            Aggregated tensor with shape ``(batch_size, n_features)``.
        """
        effective_split_size = x.shape[dim]
        if self.split_aggregation == "sum":
            # to get average
            return x.sum(dim=dim) / effective_split_size
        elif self.split_aggregation == "sum_plus":
            # to get average after softplus
            return F.softplus(x).sum(dim=dim) / effective_split_size
        elif self.split_aggregation == "logsumexp":
            # to cancel the effect of n_splits
            return torch.logsumexp(x, dim=dim) - math.log(effective_split_size)
        elif self.split_aggregation == "max":
            return torch.amax(x, dim=dim)
        else:
            raise NotImplementedError()

    def _apply_last_layer(
        self,
        last_tensor: torch.Tensor,
        cat_full_tensor: torch.Tensor | None,
        cont_full_tensor: torch.Tensor | None,
        reconstruction_indices: torch.Tensor | None = None,
    ) -> tuple[dict[str, torch.Tensor], dict[str, torch.Tensor]]:
        """Apply the last layer to the latent space.

        Parameters
        ----------
        last_tensor
            The output of the last layer.
        cat_full_tensor
            The categorical covariates.
        cont_full_tensor
            The continuous covariates.
        reconstruction_indices
            The indices of features to reconstruct.

        Returns
        -------
        tuple
            (distribution, parameters, original_parameters) where:
            - distribution: The gene expression distribution
            - parameters: Processed parameters for the distribution
            - original_parameters: Raw parameters before processing (for example pooling)
        """
        batch_size = last_tensor.shape[0]

        original_params = {}
        params = {}
        for param_name, param_info in self.gene_likelihood_module.parameters.items():
            param_net = self.params_nets[param_name]
            if param_info.startswith("fixed="):
                if reconstruction_indices is None:
                    params[param_name] = param_net.reshape(1, 1).expand(batch_size, self.n_output)
                elif reconstruction_indices.dim() == 1:
                    params[param_name] = param_net.reshape(1, 1).expand(
                        batch_size, reconstruction_indices.shape[0]
                    )
                else:
                    raise NotImplementedError()
                original_params[param_name] = params[param_name]
            elif param_info == "no_transformation":
                param_value = param_net(
                    last_tensor, cat_full_tensor, output_subset_indices=reconstruction_indices
                )
                original_params[param_name] = param_value
                if self.n_split > 1:
                    params[param_name] = self._aggregate_split(param_value, dim=-2)
                else:
                    params[param_name] = param_value
            elif param_info == "per_feature":
                param_value = param_net
                if reconstruction_indices is None:
                    pass
                elif reconstruction_indices.dim() == 1:
                    param_value = param_value[reconstruction_indices]
                else:
                    raise NotImplementedError()

                if param_value.dim() == 1:
                    param_value = param_value.unsqueeze(0).expand(batch_size, -1)
                params[param_name] = param_value
                original_params[param_name] = params[param_name]
            else:
                raise NotImplementedError()

        return params, original_params

    def forward(
        self,
        z: torch.Tensor,
        cat_full_tensor: torch.Tensor | None,
        cont_full_tensor: torch.Tensor | None,
        library: torch.Tensor,
        reconstruction_indices: torch.Tensor | None = None,
        return_original_params: bool = False,
    ) -> tuple[Any, dict[str, torch.Tensor], dict[str, torch.Tensor]]:
        """Forward computation on ``z``.

        Parameters
        ----------
        z
            Tensor of latent values with shape ``(batch_size, n_input)``.
        cat_full_tensor
            Tensor containing encoding of categorical variables of size n_batch x n_total_cat.
        cont_full_tensor
            Tensor of continuous covariate(s) for this sample.
        library
            Library size information.
        reconstruction_indices
            Indices of features to reconstruct.
        return_original_params
            Whether to return the original parameters before processing (for example pooling).

        Returns
        -------
        tuple
            (distribution, parameters, original_parameters) where:
            - distribution: The gene expression distribution
            - parameters: Processed parameters for the distribution
            - original_parameters: Raw parameters before processing (for example pooling)
        """
        # Apply split transformation
        z = self._apply_split_transformation(z)

        if cont_full_tensor is not None:
            if self.n_split > 1:
                cont_full_tensor = cont_full_tensor.unsqueeze(1).expand(-1, z.shape[1], -1)
            z = torch.cat((z, cont_full_tensor), dim=-1)

        last_tensor = z
        if self.px_shared_decoder is not None:
            last_tensor = self.px_shared_decoder(last_tensor, cat_full_tensor)
        if self.last_layer_gradient_scaler is not None:
            last_tensor = self.last_layer_gradient_scaler(last_tensor)

        params, original_params = self._apply_last_layer(
            last_tensor,
            cat_full_tensor,
            cont_full_tensor,
            reconstruction_indices=reconstruction_indices,
        )

        # Note this logic:
        px_dist = self.gene_likelihood_module.dist(parameters=params, lib_y=library)
        if not return_original_params:
            original_params = None
        return px_dist, params, original_params
