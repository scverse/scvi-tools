from __future__ import annotations

from collections import OrderedDict
from typing import Literal

import torch
from torch.distributions import Normal
from torch.nn import (
    BatchNorm1d,
    Dropout,
    LayerNorm,
    Linear,
    Module,
    Parameter,
    ReLU,
    Sequential,
)


class Embedding(Module):
    """Module for obtaining embedding of categorical covariates

    Parameters
    ----------
    size
        N categories
    cov_embed_dims
        Dimensions of embedding
    normalize
        Apply layer normalization
    """

    def __init__(self, size: int, cov_embed_dims: int = 10, normalize: bool = True):
        super().__init__()

        self.normalize = normalize

        self.embedding = torch.nn.Embedding(size, cov_embed_dims)

        if self.normalize:
            # TODO this could probably be implemented more efficiently as embed gives same result for every sample in
            #  a give class. However, if we have many balanced classes there wont be many repetitions within minibatch
            self.layer_norm = torch.nn.LayerNorm(cov_embed_dims, elementwise_affine=False)

    def forward(self, x):
        x = self.embedding(x)
        if self.normalize:
            x = self.layer_norm(x)

        return x


class EncoderDecoder(Module):
    """Module that can be used as probabilistic encoder or decoder

    Based on inputs and optional covariates predicts output mean and var

    Parameters
    ----------
    n_input
        The dimensionality of the main input
    n_output
        The dimensionality of the output
    n_cov
        Dimensionality of covariates.
        If there are no cov this should be set to None -
        in this case cov will not be used.
    n_hidden
        The number of fully-connected hidden layers
    n_layers
        Number of hidden layers
    var_mode
        How to compute variance from model outputs, see :class:`~scvi.external.sysvi.VarEncoder`
        'sample_feature' - learn per sample and feature
        'feature' - learn per feature, constant across samples
    sample
        Return samples from predicted distribution
    kwargs
        Passed to :class:`~scvi.external.sysvi.Layers`
    """

    def __init__(
        self,
        n_input: int,
        n_output: int,
        n_cov: int,
        n_hidden: int = 256,
        n_layers: int = 3,
        var_mode: Literal["sample_feature", "feature"] = "feature",
        sample: bool = False,
        **kwargs,
    ):
        super().__init__()
        self.sample = sample

        self.decoder_y = Layers(
            n_in=n_input,
            n_cov=n_cov,
            n_out=n_hidden,
            n_hidden=n_hidden,
            n_layers=n_layers,
            **kwargs,
        )

        self.mean_encoder = Linear(n_hidden, n_output)
        self.var_encoder = VarEncoder(n_hidden, n_output, mode=var_mode)

    def forward(self, x: torch.Tensor, cov: torch.Tensor | None = None):
        y = self.decoder_y(x=x, cov=cov)
        # TODO better handling of inappropriate edge-case values than nan_to_num or at least warn
        y_m = torch.nan_to_num(self.mean_encoder(y))
        y_v = self.var_encoder(y)

        outputs = {"y_m": y_m, "y_v": y_v}

        if self.sample:
            y = Normal(y_m, y_v.sqrt()).rsample()
            outputs["y"] = y

        return outputs


class Layers(Module):
    """A helper class to build fully-connected layers for a neural network.

    Adapted from scVI FCLayers to use covariates more flexibly

    Parameters
    ----------
    n_in
        The dimensionality of the main input
    n_out
        The dimensionality of the output
    n_cov
        Dimensionality of covariates.
        If there are no cov this should be set to None -
        in this case cov will not be used.
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
        Whether to inject covariates in each layer, or just the first.
    activation_fn
        Which activation function to use
    """

    def __init__(
        self,
        n_in: int,
        n_out: int,
        n_cov: int | None = None,
        n_layers: int = 1,
        n_hidden: int = 128,
        dropout_rate: float = 0.1,
        use_batch_norm: bool = True,
        use_layer_norm: bool = False,
        use_activation: bool = True,
        bias: bool = True,
        inject_covariates: bool = True,
        activation_fn: Module = ReLU,
    ):
        super().__init__()

        self.inject_covariates = inject_covariates
        self.n_cov = n_cov if n_cov is not None else 0

        layers_dim = [n_in] + (n_layers - 1) * [n_hidden] + [n_out]

        self.fc_layers = Sequential(
            OrderedDict(
                [
                    (
                        f"Layer {i}",
                        Sequential(
                            Linear(
                                n_in + self.n_cov * self.inject_into_layer(i),
                                n_out,
                                bias=bias,
                            ),
                            # non-default params come from defaults in original Tensorflow implementation
                            BatchNorm1d(n_out, momentum=0.01, eps=0.001)
                            if use_batch_norm
                            else None,
                            LayerNorm(n_out, elementwise_affine=False) if use_layer_norm else None,
                            activation_fn() if use_activation else None,
                            Dropout(p=dropout_rate) if dropout_rate > 0 else None,
                        ),
                    )
                    for i, (n_in, n_out) in enumerate(zip(layers_dim[:-1], layers_dim[1:]))
                ]
            )
        )

    def inject_into_layer(self, layer_num) -> bool:
        """Helper to determine if covariates should be injected."""
        user_cond = layer_num == 0 or (layer_num > 0 and self.inject_covariates)
        return user_cond

    def set_online_update_hooks(self, hook_first_layer=True):
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
                if isinstance(layer, Linear):
                    if self.inject_into_layer(i):
                        w = layer.weight.register_hook(_hook_fn_weight)
                    else:
                        w = layer.weight.register_hook(_hook_fn_zero_out)
                    self.hooks.append(w)
                    b = layer.bias.register_hook(_hook_fn_zero_out)
                    self.hooks.append(b)

    def forward(self, x: torch.Tensor, cov: torch.Tensor | None = None):
        """
        Forward computation on ``x``.

        Parameters
        ----------
        x
            tensor of values with shape ``(n_in,)``
        cov
            tensor of covariate values with shape ``(n_cov,)`` or None

        Returns
        -------
        py:class:`torch.Tensor`
            tensor of shape ``(n_out,)``

        """
        for i, layers in enumerate(self.fc_layers):
            for layer in layers:
                if layer is not None:
                    if isinstance(layer, BatchNorm1d):
                        if x.dim() == 3:
                            x = torch.cat([(layer(slice_x)).unsqueeze(0) for slice_x in x], dim=0)
                        else:
                            x = layer(x)
                    else:
                        # Injection of covariates
                        if (
                            self.n_cov > 0
                            and isinstance(layer, Linear)
                            and self.inject_into_layer(i)
                        ):
                            x = torch.cat((x, cov), dim=-1)
                        x = layer(x)
        return x


class VarEncoder(Module):
    """Encode variance (strictly positive).

    Parameters
    ----------
    n_input
        Number of input dimensions, used if mode is sample_feature
    n_output
        Number of variances to predict
    mode
        How to compute var
        'sample_feature' - learn per sample and feature
        'feature' - learn per feature, constant across samples
    """

    def __init__(
        self,
        n_input: int,
        n_output: int,
        mode: Literal["sample_feature", "feature", "linear"],
    ):
        super().__init__()

        self.eps = 1e-4
        self.mode = mode
        if self.mode == "sample_feature":
            self.encoder = Linear(n_input, n_output)
        elif self.mode == "feature":
            self.var_param = Parameter(torch.zeros(1, n_output))
        else:
            raise ValueError("Mode not recognised.")
        self.activation = torch.exp

    def forward(self, x: torch.Tensor):
        """Forward pass through model

        Parameters
        ----------
        x
            Used to encode var if mode is sample_feature; dim = n_samples x n_input

        Returns
        -------
        Predicted var
        """
        # Force to be non nan - TODO come up with better way to do so
        if self.mode == "sample_feature":
            v = self.encoder(x)
            v = (self.activation(v) + self.eps)  # Ensure that var is strictly positive
        elif self.mode == "feature":
            v = self.var_param.expand(x.shape[0], -1)  # Broadcast to input size
            v = (self.activation(v) + self.eps)  # Ensure that var is strictly positive
        return v
