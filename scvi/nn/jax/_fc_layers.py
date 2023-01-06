from typing import List, Literal, Optional

import jax.numpy as jnp
from flax import linen as nn
from flax.linen.initializers import variance_scaling

BATCH_NORM_KWARGS = {
    "momentum": 0.99,
    "epsilon": 0.001,
}
LAYER_NORM_KWARGS = {
    "use_bias": False,  # replicates pytorch elementwise_affine=False
    "use_scale": False,  # replicates pytorch elementwise_affine=False
}


class Dense(nn.Dense):
    """Jax dense layer with defaults to replicate PyTorch initialization."""

    def __init__(self, *args, **kwargs):
        scale = 1 / 3
        kernel_init = variance_scaling(scale, "fan_in", "uniform")
        # bias init can't see input shape so don't include here
        kwargs.update({"kernel_init": kernel_init})
        super().__init__(*args, **kwargs)


class JaxFCLayers(nn.Module):
    """Fully-connected layers with a Jax backend."""

    n_input: int
    layers_dim: List[int]
    bias: bool
    dropout_rate: float
    norm: Literal["batch", "layer", "group", None]
    norm_kwargs: dict
    activation: str
    activation_kwargs: dict
    residual: bool
    n_cat_list: List[int]
    inject_covariates: bool
    training: Optional[bool]

    def setup(self):
        """Setup the module."""
        self._norm = None
        if self.norm == "batch":
            self._norm = nn.BatchNorm
            self._norm_kwargs = BATCH_NORM_KWARGS.copy()
            self._norm_kwargs.update(self.norm_kwargs)
        elif self.norm == "layer":
            self._norm = nn.LayerNorm
            self._norm_kwargs = LAYER_NORM_KWARGS.copy()
            self._norm_kwargs.update(self.norm_kwargs)
        elif self.norm == "group":
            self._norm = nn.GroupNorm

    def _forward_block(
        self,
        x: jnp.ndarray,
        block_num: int,
        cat_list: Optional[List[int]] = None,
        is_eval: bool = None,
    ) -> jnp.ndarray:
        """Forward pass through a single fully-connected block."""
        # TODO: inject covariates
        x = Dense(self.layers_dim[block_num], use_bias=self.bias)(x)
        x = (
            self._norm(**self._norm_kwargs)(x, deterministic=is_eval)
            if self.norm
            else x
        )
        x = self._activation(x, **self._activation_kwargs) if self.activation else x
        x = (
            nn.Dropout(rate=self.dropout_rate)(x, deterministic=is_eval)
            if self.dropout_rate > 0
            else x
        )
        return x

    @nn.compact
    def __call__(
        self,
        x: jnp.ndarray,
        cat_list: Optional[List[int]] = None,
        training: Optional[bool] = None,
    ) -> jnp.ndarray:
        """
        Forward pass through the module.

        Parameters
        ----------
        x
            Input tensor.
        cat_list
            List of categorical membership(s).

        """
        training = nn.merge_param("training", self.training, training)
        is_eval = not training

        x = self._forward_block(x, 0, cat_list, is_eval)

        for block_num in range(1, len(self.layers_dim)):
            x_ = self._forward_block(x, block_num, is_eval=is_eval)
            x = x_ + x if self.residual else x_

        return x

    def _inject_into_layer(self, layer_num) -> bool:
        """Determines whether covariates should be injected into a given layer."""
        return layer_num == 0 or self.inject_covariates
