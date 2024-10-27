from __future__ import annotations

import collections
import warnings
from collections.abc import Iterable
from typing import Literal, Callable

import numpy as np
import torch
from torch import nn
from torch.distributions import Normal
from torch.nn import (
    Linear,
    Module,
    Parameter,
)


class EncoderDecoder(Module):
    """Module that can be used as probabilistic encoder or decoder.

    Based on inputs and optional covariates predicts output mean and variance.

    Parameters
    ----------
    n_input
        The dimensionality of the main input.
    n_output
        The dimensionality of the output.
    n_cat_list
        A list containing the number of categories for each covariate.
    n_cont
        The dimensionality of the continuous covariates.
    n_hidden
        The number of nodes per hidden layer.
    n_layers
        The number of hidden layers.
    var_mode
        How to compute variance from model outputs, see :class:`~scvi.external.sysvi.VarEncoder`.
        One of the following:
        * ```'sample_feature'``` - learn variance per sample and feature.
        * ```'feature'``` - learn variance per feature, constant across samples.
    var_activation
        Function used to ensure positivity of the variance.
        Defaults to :meth:`torch.exp`.
    sample
        Return samples from predicted distribution.
    kwargs
         Passed to :class:`~scvi.external.sysvi.Layers`.
    """

    def __init__(
        self,
        n_input: int,
        n_output: int,
        n_cat_list: list[int],
        n_cont: int,
        n_hidden: int = 256,
        n_layers: int = 3,
        var_mode: Literal["sample_feature", "feature"] = "feature",
        var_activation: Callable | None = None,
        sample: bool = False,
        **kwargs,
    ):
        super().__init__()
        self.sample = sample

        self.decoder_y = FCLayers(
            n_in=n_input,
            n_cat_list=n_cat_list,
            n_cont=n_cont,
            n_out=n_hidden,
            n_hidden=n_hidden,
            n_layers=n_layers,
            **kwargs,
        )

        self.mean_encoder = Linear(n_hidden, n_output)
        self.var_encoder = VarEncoder(n_hidden, n_output, mode=var_mode, activation=var_activation)

    def forward(
        self,
        x: torch.Tensor,
        cont: torch.Tensor | None = None,
        cat_list: list[torch.Tensor] | None = None,
    ) -> dict[str, torch.Tensor]:
        """Forward pass.

        Parameters
        ----------
        x
            Main input (i.e. expression for encoder or latent embedding for decoder.).
            dim = n_samples * n_input
        cont
            Optional continuous covariates.
            dim = n_samples * n_cont
        cat_list
            List of optional categorical covariates.
            Will be one hot encoded in `~scvi.nn.FCLayers`.
            Each list entry is of dim = n_samples * 1

        Returns
        -------
        Predicted mean (``'y_m'``) and variance (``'y_v'``) and
        optionally samples (``'y'``) from normal distribution
        parametrized with the predicted parameters.
        """
        y = self.decoder_y(x=x, cont=cont, cat_list=cat_list)
        y_m = self.mean_encoder(y)
        if y_m.isnan().any() or y_m.isinf().any():
            warnings.warn("Predicted mean contains nan or inf values. Setting to numerical.")
            y_m = torch.nan_to_num(y_m)
        y_v = self.var_encoder(y)

        outputs = {"y_m": y_m, "y_v": y_v}

        if self.sample:
            y = Normal(y_m, y_v.sqrt()).rsample()
            outputs["y"] = y

        return outputs


class FCLayers(nn.Module):
    """A helper class to build fully-connected layers for a neural network.

    FCLayers class of scvi-tools adapted to also inject continous covariates.

    The only adaptation is addition of `n_cont` parameter in init and `cont` in forward,
    with the associated handling of the two.
    The forward method signature is changed to account for optional `cont`.

    Parameters
    ----------
    n_in
        The dimensionality of the input
    n_out
        The dimensionality of the output
    n_cat_list
        The number of categorical covariates and
        the number of category levels.
        A list containing, for each covariate of interest,
        the number of categories. Each covariate will be
        included using a one-hot encoding.
    n_cont
        The number of continuous covariates.
    n_layers
        The number of fully-connected hidden layers
    n_hidden
        The number of nodes per hidden layer
    dropout_rate
        Dropout rate to apply to each of the hidden layers
    use_batch_norm
        Whether to have `BatchNorm` layers or not
    use_layer_norm
        Whether to have `LayerNorm` layers or not
    use_activation
        Whether to have layer activation or not
    bias
        Whether to learn bias in linear layers or not
    inject_covariates
        Whether to inject covariates in each layer, or just the first (default).
    activation_fn
        Which activation function to use
    """

    def __init__(
        self,
        n_in: int,
        n_out: int,
        n_cat_list: Iterable[int] = None,
        n_cont: int = 0,
        n_layers: int = 1,
        n_hidden: int = 128,
        dropout_rate: float = 0.1,
        use_batch_norm: bool = True,
        use_layer_norm: bool = False,
        use_activation: bool = True,
        bias: bool = True,
        inject_covariates: bool = True,
        activation_fn: nn.Module = nn.ReLU,
    ):
        super().__init__()
        self.inject_covariates = inject_covariates
        layers_dim = [n_in] + (n_layers - 1) * [n_hidden] + [n_out]

        if n_cat_list is not None:
            # n_cat = 1 will be ignored
            self.n_cat_list = [n_cat if n_cat > 1 else 0 for n_cat in n_cat_list]
        else:
            self.n_cat_list = []

        self.n_cov = sum(self.n_cat_list) + n_cont
        self.fc_layers = nn.Sequential(
            collections.OrderedDict(
                [
                    (
                        f"Layer {i}",
                        nn.Sequential(
                            nn.Linear(
                                n_in + self.n_cov * self.inject_into_layer(i),
                                n_out,
                                bias=bias,
                            ),
                            # non-default params come from defaults in original Tensorflow
                            # implementation
                            nn.BatchNorm1d(n_out, momentum=0.01, eps=0.001)
                            if use_batch_norm
                            else None,
                            nn.LayerNorm(n_out, elementwise_affine=False)
                            if use_layer_norm
                            else None,
                            activation_fn() if use_activation else None,
                            nn.Dropout(p=dropout_rate) if dropout_rate > 0 else None,
                        ),
                    )
                    for i, (n_in, n_out) in enumerate(
                    zip(layers_dim[:-1], layers_dim[1:], strict=True)
                )
                ]
            )
        )

    def inject_into_layer(self, layer_num) -> bool:
        """Helper to determine if covariates should be injected."""
        user_cond = layer_num == 0 or (layer_num > 0 and self.inject_covariates)
        return user_cond

    def set_online_update_hooks(self, hook_first_layer=True):
        """Set online update hooks."""
        self.hooks = []

        def _hook_fn_weight(grad):
            new_grad = torch.zeros_like(grad)
            if self.n_cov > 0:
                new_grad[:, -self.n_cov:] = grad[:, -self.n_cov:]
            return new_grad

        def _hook_fn_zero_out(grad):
            return grad * 0

        for i, layers in enumerate(self.fc_layers):
            for layer in layers:
                if i == 0 and not hook_first_layer:
                    continue
                if isinstance(layer, nn.Linear):
                    if self.inject_into_layer(i):
                        w = layer.weight.register_hook(_hook_fn_weight)
                    else:
                        w = layer.weight.register_hook(_hook_fn_zero_out)
                    self.hooks.append(w)
                    b = layer.bias.register_hook(_hook_fn_zero_out)
                    self.hooks.append(b)

    def forward(
        self, x: torch.Tensor, cont: torch.Tensor | None = None, cat_list: list | None = None
    ) -> torch.Tensor:
        """Forward computation on ``x``.

        Parameters
        ----------
        x
            tensor of values with shape ``(n_in,)``
        cont
            continuous covariates for this sample,
            tensor of values with shape ``(n_cont,)``
        cat_list
            list of category membership(s) for this sample

        Returns
        -------
        :class:`torch.Tensor`
            tensor of shape ``(n_out,)``
        """
        one_hot_cat_list = []  # for generality in this list many indices useless.
        cont_list = [cont] if cont is not None else []
        cat_list = cat_list or []

        if len(self.n_cat_list) > len(cat_list):
            raise ValueError("nb. categorical args provided doesn't match init. params.")
        for n_cat, cat in zip(self.n_cat_list, cat_list, strict=False):
            if n_cat and cat is None:
                raise ValueError("cat not provided while n_cat != 0 in init. params.")
            if n_cat > 1:  # n_cat = 1 will be ignored - no additional information
                if cat.size(1) != n_cat:
                    one_hot_cat = nn.functional.one_hot(cat.squeeze(-1), n_cat)
                else:
                    one_hot_cat = cat  # cat has already been one_hot encoded
                one_hot_cat_list += [one_hot_cat]
        for i, layers in enumerate(self.fc_layers):
            for layer in layers:
                if layer is not None:
                    if isinstance(layer, nn.BatchNorm1d):
                        if x.dim() == 3:
                            x = torch.cat([(layer(slice_x)).unsqueeze(0) for slice_x in x], dim=0)
                        else:
                            x = layer(x)
                    else:
                        if isinstance(layer, nn.Linear) and self.inject_into_layer(i):
                            if x.dim() == 3:
                                cov_list_layer = [
                                    o.unsqueeze(0).expand((x.size(0), o.size(0), o.size(1)))
                                    for o in one_hot_cat_list
                                ]
                            else:
                                cov_list_layer = one_hot_cat_list
                            x = torch.cat((x, *cov_list_layer, *cont_list), dim=-1)
                        x = layer(x)
        return x


class VarEncoder(Module):
    """Encode variance (strictly positive).

    Parameters
    ----------
    n_input
        Number of input dimensions.
        Used if mode is ``'sample_feature'`` to construct a network predicting
        variance from input features.
    n_output
        Number of variances to predict, matching the desired number of features
        (e.g. latent dimensions for variational encoding or output features
        for variational decoding).
    mode
        How to compute variance.
        One of the following:
        * ``'sample_feature'`` - learn variance per sample and feature.
        * ``'feature'`` - learn variance per feature, constant across samples.
    """

    def __init__(
        self,
        n_input: int,
        n_output: int,
        mode: Literal["sample_feature", "feature", "linear"],
        activation: Callable | None = None,
    ):
        super().__init__()

        self.clip_exp = np.log(torch.finfo(torch.get_default_dtype()).max) - 1e-4
        self.mode = mode
        if self.mode == "sample_feature":
            self.encoder = Linear(n_input, n_output)
        elif self.mode == "feature":
            self.var_param = Parameter(torch.zeros(1, n_output))
        else:
            raise ValueError("Mode not recognised.")
        self.activation = torch.exp if activation is None else activation

    def forward(
        self,
        x: torch.Tensor,
    ) -> torch.Tensor:
        """Forward pass through model.

        Parameters
        ----------
        x
            Used to encode variance if mode is ``'sample_feature'``.
            dim = n_samples x n_input

        Returns
        -------
        Predicted variance
        dim = n_samples * 1
        """
        if self.mode == "sample_feature":
            v = self.encoder(x)
        elif self.mode == "feature":
            v = self.var_param.expand(x.shape[0], -1)  # Broadcast to input size

        # Ensure that var is strictly positive via exp - Bring back to non-log scale
        # Clip to range that will not be inf after exp
        v = torch.clip(v, min=-self.clip_exp, max=self.clip_exp)
        v = self.activation(v)
        if v.isnan().any():
            warnings.warn("Predicted variance contains nan values. Setting to 0.")
            v = torch.nan_to_num(v)

        return v
