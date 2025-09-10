from __future__ import annotations

from typing import TYPE_CHECKING

import flax.linen as nn
import jax
import jax.numpy as jnp
import numpyro.distributions as dist

if TYPE_CHECKING:
    from collections.abc import Callable
    from typing import Any, Literal

PYTORCH_DEFAULT_SCALE = 1 / 3


class Dense(nn.DenseGeneral):
    """Dense layer.

    Uses a custom initializer for the kernel to replicate the default PyTorch behavior.
    """

    def __init__(self, *args, **kwargs):
        from flax.linen.initializers import variance_scaling

        _kwargs = {"kernel_init": variance_scaling(PYTORCH_DEFAULT_SCALE, "fan_in", "uniform")}
        _kwargs.update(kwargs)

        super().__init__(*args, **_kwargs)


class ResnetBlock(nn.Module):
    """Resnet block.

    Consists of the following operations:

    1. :class:`~flax.linen.Dense`
    2. :class:`~flax.linen.LayerNorm`
    3. Activation function specified by ``internal_activation``
    4. Skip connection if ``n_in`` is equal to ``n_hidden``, otherwise a :class:`~flax.linen.Dense`
        layer is applied to the input before the skip connection to match features.
    5. :class:`~flax.linen.Dense`
    6. :class:`~flax.linen.LayerNorm`
    7. Activation function specified by ``output_activation``

    Parameters
    ----------
    n_out
        Number of output units.
    n_hidden
        Number of hidden units.
    internal_activation
        Activation function to use after the first :class:`~flax.linen.Dense` layer.
    output_activation
        Activation function to use after the last :class:`~flax.linen.Dense` layer.
    training
        Whether the model is in training mode.
    """

    n_out: int
    n_hidden: int = 128
    internal_activation: Callable[[jax.typing.ArrayLike], jax.Array] = nn.relu
    output_activation: Callable[[jax.typing.ArrayLike], jax.Array] = nn.relu
    training: bool | None = None

    @nn.compact
    def __call__(self, inputs: jax.typing.ArrayLike, training: bool | None = None) -> jax.Array:
        training = nn.merge_param("training", self.training, training)
        h = Dense(self.n_hidden)(inputs)
        h = nn.LayerNorm()(h)
        h = self.internal_activation(h)
        if inputs.shape[-1] != self.n_hidden:
            h = h + Dense(self.n_hidden)(inputs)
        else:
            h = h + inputs
        h = Dense(self.n_out)(h)
        h = nn.LayerNorm()(h)
        return self.output_activation(h)


class MLP(nn.Module):
    """Multi-layer perceptron with resnet blocks.

    Applies ``n_layers`` :class:`~ResnetBlock` blocks to the input, followed by a
    :class:`~flax.linen.Dense` layer to project to the output dimension.

    Parameters
    ----------
    n_out
        Number of output units.
    n_hidden
        Number of hidden units.
    n_layers
        Number of resnet blocks.
    activation
        Activation function to use.
    training
        Whether the model is in training mode.
    """

    n_out: int
    n_hidden: int = 128
    n_layers: int = 1
    activation: Callable[[jax.typing.ArrayLike], jax.Array] = nn.relu
    training: bool | None = None

    @nn.compact
    def __call__(self, inputs: jax.typing.ArrayLike, training: bool | None = None) -> jax.Array:
        training = nn.merge_param("training", self.training, training)
        h = inputs
        for _ in range(self.n_layers):
            h = ResnetBlock(
                n_out=self.n_hidden,
                internal_activation=self.activation,
                output_activation=self.activation,
            )(h, training=training)
        return Dense(self.n_out)(h)


class NormalDistOutputNN(nn.Module):
    """Fully-connected neural net parameterizing a normal distribution.

    Applies ``n_layers`` :class:`~ResnetBlock` blocks to the input, followed by a
    :class:`~flax.linen.Dense` layer for the mean and a :class:`~flax.linen.Dense` and
    :func:`~flax.linen.softplus` layer for the scale.

    Parameters
    ----------
    n_out
        Number of output units.
    n_hidden
        Number of hidden units.
    n_layers
        Number of resnet blocks.
    scale_eps
        Numerical stability constant added to the scale of the normal distribution.
    """

    n_out: int
    n_hidden: int = 128
    n_layers: int = 1
    scale_eps: float = 1e-5
    training: bool | None = None

    @nn.compact
    def __call__(self, inputs: jax.typing.ArrayLike, training: bool | None = None) -> dist.Normal:
        training = nn.merge_param("training", self.training, training)
        h = inputs
        for _ in range(self.n_layers):
            h = ResnetBlock(n_out=self.n_hidden)(h, training=training)
        mean = Dense(self.n_out)(h)
        scale = nn.Sequential([Dense(self.n_out), nn.softplus])(h)
        return dist.Normal(mean, scale + self.scale_eps)


class ConditionalNormalization(nn.Module):
    """Condition-specific normalization.

    Applies either batch normalization or layer normalization to the input, followed by
    condition-specific scaling (``gamma``) and shifting (``beta``).

    Parameters
    ----------
    n_features
        Number of features.
    n_conditions
        Number of conditions.
    training
        Whether the model is in training mode.
    normalization_type
        Type of normalization to apply. Must be one of ``"batch", "layer"``.
    """

    n_features: int
    n_conditions: int
    training: bool | None = None
    normalization_type: Literal["batch", "layer"] = "layer"

    @staticmethod
    def _gamma_initializer() -> jax.nn.initializers.Initializer:
        def init(key: jax.random.KeyArray, shape: tuple, dtype: Any = jnp.float_) -> jax.Array:
            weights = jax.random.normal(key, shape, dtype) * 0.02 + 1
            return weights

        return init

    @staticmethod
    def _beta_initializer() -> jax.nn.initializers.Initializer:
        def init(key: jax.random.KeyArray, shape: tuple, dtype: Any = jnp.float_) -> jax.Array:
            del key
            weights = jnp.zeros(shape, dtype=dtype)
            return weights

        return init

    @nn.compact
    def __call__(
        self,
        x: jax.typing.ArrayLike,
        condition: jax.typing.ArrayLike,
        training: bool | None = None,
    ) -> jax.Array:
        training = nn.merge_param("training", self.training, training)

        if self.normalization_type == "batch":
            x = nn.BatchNorm(use_bias=False, use_scale=False)(x, use_running_average=not training)
        elif self.normalization_type == "layer":
            x = nn.LayerNorm(use_bias=False, use_scale=False)(x)
        else:
            raise ValueError("`normalization_type` must be one of ['batch', 'layer'].")

        cond_int = condition.squeeze(-1).astype(int)
        gamma = nn.Embed(
            self.n_conditions,
            self.n_features,
            embedding_init=self._gamma_initializer(),
            name="gamma_conditional",
        )(cond_int)
        beta = nn.Embed(
            self.n_conditions,
            self.n_features,
            embedding_init=self._beta_initializer(),
            name="beta_conditional",
        )(cond_int)

        return gamma * x + beta


class AttentionBlock(nn.Module):
    """Attention block consisting of multi-head self-attention and MLP.

    Parameters
    ----------
    query_dim
        Dimension of the query input.
    out_dim
        Dimension of the output.
    outerprod_dim
        Dimension of the outer product.
    n_channels
        Number of channels.
    n_heads
        Number of heads.
    dropout_rate
        Dropout rate.
    n_hidden_mlp
        Number of hidden units in the MLP.
    n_layers_mlp
        Number of layers in the MLP.
    training
        Whether the model is in training mode.
    stop_gradients_mlp
        Whether to stop gradients through the MLP.
    activation
        Activation function to use.
    """

    query_dim: int
    out_dim: int
    outerprod_dim: int = 16
    n_channels: int = 4
    n_heads: int = 2
    dropout_rate: float = 0.0
    n_hidden_mlp: int = 32
    n_layers_mlp: int = 1
    training: bool | None = None
    stop_gradients_mlp: bool = False
    activation: Callable[[jax.Array], jax.Array] = nn.gelu

    @nn.compact
    def __call__(
        self,
        query_embed: jax.typing.ArrayLike,
        kv_embed: jax.typing.ArrayLike,
        training: bool | None = None,
    ) -> jax.Array:
        training = nn.merge_param("training", self.training, training)
        has_mc_samples = query_embed.ndim == 3

        query_embed_stop = (
            query_embed if not self.stop_gradients_mlp else jax.lax.stop_gradient(query_embed)
        )
        query_for_att = nn.DenseGeneral((self.outerprod_dim, 1), use_bias=False)(query_embed_stop)
        kv_for_att = nn.DenseGeneral((self.outerprod_dim, 1), use_bias=False)(kv_embed)
        eps = nn.MultiHeadDotProductAttention(
            num_heads=self.n_heads,
            qkv_features=self.n_channels * self.n_heads,
            out_features=self.n_channels,
            dropout_rate=self.dropout_rate,
            use_bias=True,
        )(inputs_q=query_for_att, inputs_kv=kv_for_att, deterministic=not training)

        if not has_mc_samples:
            eps = jnp.reshape(eps, (eps.shape[0], eps.shape[1] * eps.shape[2]))
        else:
            eps = jnp.reshape(eps, (eps.shape[0], eps.shape[1], eps.shape[2] * eps.shape[3]))

        eps_ = MLP(
            n_out=self.outerprod_dim,
            n_hidden=self.n_hidden_mlp,
            training=training,
            activation=self.activation,
        )(inputs=eps)
        inputs = jnp.concatenate([query_embed, eps_], axis=-1)
        residual = MLP(
            n_out=self.out_dim,
            n_hidden=self.n_hidden_mlp,
            n_layers=self.n_layers_mlp,
            training=training,
            activation=self.activation,
        )(inputs=inputs)
        return residual
