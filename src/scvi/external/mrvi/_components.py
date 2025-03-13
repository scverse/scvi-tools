from __future__ import annotations

from typing import TYPE_CHECKING

# import flax.linen as nn
import jax
import jax.numpy as jnp

if TYPE_CHECKING:
    from collections.abc import Callable
    from typing import Literal

    import torch

from torch import nn
from torch.distributions import Normal

PYTORCH_DEFAULT_SCALE = 1 / 3


class Dense(nn.Linear):
    def __init__(self, *args, **kwargs):
        # TODO: Might need to change the kernel initialization.
        # Not sure if default torch behavior is correct
        # otherwise can just get rid of this class and use nn.Linear itself

        super().__init__(*args, **kwargs)


class ResnetBlock(nn.Module):
    def __init__(
        self,
        # TODO: should I keep n_in or is there a functional way to do this like in flax?
        n_in: int,
        n_out: int,
        n_hidden: int = 128,
        internal_activation: Callable[[torch.Tensor], torch.Tensor] = nn.ReLU,
        output_activation: Callable[[torch.Tensor], torch.Tensor] = nn.ReLU,
        training: bool | None = None,
    ):
        super.__init__()
        self.n_in = n_in
        self.n_out = n_out
        self.n_hidden = n_hidden
        self.internal_activation = internal_activation
        self.output_activation = output_activation
        self.training = training

        # TODO: figure out what below does (taken from jax code)
        # training = nn.merge_param("training", self.training, training)

        # dense layer
        self.fc1 = Dense(in_features=n_in, out_features=n_out)

        # layer norm
        self.layer_norm1 = nn.LayerNorm(n_out)

        # internal activation

        # skip connection if n_in equal to n_hidden,
        # otherwise dense layer applied before skip connection to match features
        self.fc_match = Dense(in_features=n_in, out_features=n_hidden)

        # dense layer
        self.fc2 = Dense(in_features=n_hidden, out_features=n_out)

        # layer norm
        self.layer_norm2 = nn.LayerNorm(n_out)

        # output activation

    def forward(self, inputs: torch.Tensor) -> torch.Tensor:
        h = self.fc1(inputs)
        h = self.layer_norm1(h)
        h = self.internal_activation(h)

        if self.n_in != self.n_hidden:
            h = h + self.fc_match(inputs)
        else:
            h = h + inputs

        h = self.fc2(h)
        h = self.layer_norm2(h)
        return self.output_activation(h)


class MLP(nn.Module):
    def __init__(
        self,
        n_in: int,
        n_out: int,
        n_hidden: int = 128,
        n_layers: int = 1,
        activation: Callable[[torch.Tensor], torch.Tensor] = nn.ReLU,
        training: bool | None = None,
    ):
        super.__init__()
        self.n_in = n_in
        self.n_out = n_out
        self.n_hidden = n_hidden
        self.n_layers = n_layers
        self.activation = activation
        self.training = training

        # sequence of n_layers resnet blocks
        self.resnet_blocks = nn.Sequential(
            *[
                ResnetBlock(
                    n_in=n_in,
                    n_out=n_hidden,
                    internal_activation=activation,
                    output_activation=activation,
                )
                for _ in range(n_layers)
            ]
        )

        # dense layer to project to the output dimension
        self.fc = Dense(in_features=n_hidden, out_features=n_out)

    def forward(self, inputs: torch.Tensor, training: bool | None = None) -> torch.Tensor:
        # TODO: figure out what below is
        # training = nn.merge_param("training", self.training, training)

        h = self.resnet_blocks(inputs)
        return self.fc(h)


class NormalDistOutputNN(nn.Module):
    def __init__(
        self,
        n_in: int,
        n_out: int,
        n_hidden: int = 128,
        n_layers: int = 1,
        scale_eps: float = 1e-5,
        training: bool | None = None,
    ):
        super().__init__()
        self.n_in = n_in
        self.n_out = n_out
        self.n_hidden = n_hidden
        self.n_layers = n_layers
        self.scale_eps = scale_eps
        self.training = training

        self.resnet_blocks = nn.Sequential(
            *[
                ResnetBlock(
                    n_in=n_in,
                    n_out=n_hidden,
                )
                for _ in range(n_layers)
            ]
        )

        self.fc_mean = Dense(in_features=n_hidden, out_features=n_out)
        self.fc_scale = nn.Sequential(
            Dense(in_features=n_hidden, out_features=n_out),
            nn.Softplus(),
        )

    def forward(self, inputs: torch.Tensor, training: bool | None = None) -> Normal:
        # what is below doing?
        # training = nn.merge_param("training", self.training, training)

        h = self.resnet_blocks(inputs, training)
        mean = self.fc_mean(h)
        scale = self.fc_scale(h)
        return Normal(mean, scale + self.scale_eps)


# TODO: implement ConditionalNormalization in torch
class ConditionalNormalization(nn.Module):
    def __init__(
        self,
        n_features: int,
        n_conditions: int,
        training: bool | None = None,
        normalization_type: Literal["batch", "layer"] = "layer",
    ):
        super().__init__()
        self.n_features = n_features
        self.n_conditions = n_conditions
        self.training = training
        self.normalization_type = normalization_type

    def forward(self, x: torch.Tensor, condition: torch.Tensor, training: bool | None = None):
        # training = nn.merge_param("training", self.training, training)
        if self.normalization_type == "batch":
            x = nn.functional.batch_norm(
                x, running_mean=not training, running_var=not training, training=training
            )
        elif self.normalization_type == "layer":
            x = nn.functional.layer_norm(x, training=training)
        else:
            raise ValueError("`normalization_type` must be one of ['batch', 'layer'].")

        # TODO: figure out how to initialize embeddings in torch
        cond_int = condition.squeeze(-1).astype(int)
        gamma = nn.Embedding(self.n_conditions, self.n_features)(cond_int)
        beta = nn.Embedding(self.n_conditions, self.n_features)(cond_int)

        return gamma * x + beta


# TODO: implement AttentionBlock in torch


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
