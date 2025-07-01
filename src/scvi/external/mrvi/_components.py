from __future__ import annotations

from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from collections.abc import Callable
    from typing import Literal

import torch
from torch import nn
from torch.distributions import Normal

PYTORCH_DEFAULT_SCALE = 1 / 3


class Dense(nn.Linear):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)


class ResnetBlock(nn.Module):
    def __init__(
        self,
        # TODO: should I keep n_in or is there a functional way to do this like in flax?
        n_in: int,
        n_out: int,
        n_hidden: int = 128,
        internal_activation: Callable[[torch.Tensor], torch.Tensor] = nn.functional.relu,
        output_activation: Callable[[torch.Tensor], torch.Tensor] = nn.functional.relu,
        training: bool | None = None,
    ):
        super().__init__()
        self.n_in = n_in
        self.n_out = n_out
        self.n_hidden = n_hidden
        self.internal_activation = internal_activation
        self.output_activation = output_activation
        self.training = training

        # dense layer
        self.fc1 = Dense(in_features=n_in, out_features=n_hidden)

        # layer norm
        self.layer_norm1 = nn.LayerNorm(n_hidden)

        # internal activation

        # skip connection if n_in equal to n_hidden,
        # otherwise dense layer applied before skip connection to match features
        self.fc_match = Dense(in_features=n_in, out_features=n_hidden)

        # dense layer
        self.fc2 = Dense(in_features=n_hidden, out_features=n_out)

        # layer norm
        self.layer_norm2 = nn.LayerNorm(n_out)

    def forward(self, inputs: torch.Tensor, training: bool | None = None) -> torch.Tensor:
        self.training = training  # is this sufficients compared to jax version?
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
        super().__init__()
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
        self.training = training  # is this sufficient compared to jax version?

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

        self.resnet_blocks = nn.ModuleList(
            [
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
        h = inputs
        for block in self.resnet_blocks:
            h = block(h, training)
        mean = self.fc_mean(h)
        scale = self.fc_scale(h)
        return Normal(mean, scale + self.scale_eps)


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

        # Initialize embedding layers in __init__ instead of forward
        self.gamma_embedding = nn.Embedding(self.n_conditions, self.n_features)
        self.beta_embedding = nn.Embedding(self.n_conditions, self.n_features)

        # Initialize normalization layers
        if self.normalization_type == "batch":
            self.norm_layer = nn.BatchNorm1d(
                self.n_features, affine=False, track_running_stats=True
            )
        elif self.normalization_type == "layer":
            self.norm_layer = nn.LayerNorm(self.n_features, elementwise_affine=False)
        else:
            raise ValueError("`normalization_type` must be one of ['batch', 'layer'].")

    def forward(self, x: torch.Tensor, condition: torch.Tensor, training: bool | None = None):
        # Use pre-initialized normalization layer
        if self.normalization_type == "batch":
            # For BatchNorm, we need to set training mode
            if training is not None:
                self.norm_layer.training = training
            x = self.norm_layer(x)
        else:  # layer norm
            x = self.norm_layer(x)

        # Use pre-initialized embedding layers
        cond_int = condition.squeeze(-1).to(torch.int64)  # TODO: check if 64 or 32
        gamma = self.gamma_embedding(cond_int)
        beta = self.beta_embedding(cond_int)

        return gamma * x + beta


class AttentionBlock(nn.Module):
    def __init__(
        self,
        query_dim: int,
        kv_dim: int,  # TODO: added this for torch version, could be wrong
        out_dim: int,
        outerprod_dim: int = 16,
        n_channels: int = 4,
        n_heads: int = 2,
        dropout_rate: float = 0.0,
        n_hidden_mlp: int = 32,
        n_layers_mlp: int = 1,
        training: bool | None = None,
        stop_gradients_mlp: bool = False,
        activation: Callable[[torch.Tensor], torch.Tensor] = nn.functional.gelu,
    ):
        super().__init__()
        self.query_dim = query_dim
        self.kv_dim = kv_dim
        self.out_dim = out_dim
        self.outerprod_dim = outerprod_dim
        self.n_channels = n_channels
        self.n_heads = n_heads
        self.dropout_rate = dropout_rate
        self.n_hidden_mlp = n_hidden_mlp
        self.n_layers_mlp = n_layers_mlp
        self.training = training
        self.stop_gradients_mlp = stop_gradients_mlp
        self.activation = activation

        self.query_proj = nn.Linear(in_features=query_dim, out_features=outerprod_dim, bias=False)
        self.embed_dim_proj_query = nn.Linear(
            in_features=1, out_features=n_channels * n_heads, bias=False
        )
        self.embed_dim_proj_kv = nn.Linear(
            in_features=1, out_features=n_channels * n_heads, bias=False
        )
        self.kv_proj = nn.Linear(in_features=kv_dim, out_features=outerprod_dim, bias=False)
        self.attention = nn.MultiheadAttention(
            embed_dim=n_channels * n_heads,
            num_heads=n_heads,
            dropout=dropout_rate,
            batch_first=True,
        )
        self.mlp_eps = MLP(
            n_in=self.outerprod_dim * n_channels * n_heads,
            n_out=outerprod_dim,
            n_hidden=n_hidden_mlp,
            n_layers=n_layers_mlp,
            training=training,
            activation=activation,
        )
        self.mlp_residual = MLP(
            n_in=outerprod_dim + query_dim,
            n_out=out_dim,
            n_hidden=n_hidden_mlp,
            n_layers=n_layers_mlp,
            training=training,
            activation=activation,
        )

    def forward(
        self, query_embed: torch.Tensor, kv_embed: torch.Tensor, training: bool | None = None
    ) -> torch.Tensor:
        self.training = training
        has_mc_samples = query_embed.ndim == 3

        if self.stop_gradients_mlp:
            query_embed_stop = query_embed.detach()
        else:
            query_embed_stop = query_embed  # (batch_size, query_dim)

        # Below, the second projection was not needed in the original JAX code,
        # but it is needed in the PyTorch version to match the dimensions
        # of the query and key-value embeddings for the attention mechanism.
        query_for_att = self.query_proj(query_embed_stop).unsqueeze(
            -1
        )  # (batch_size, outerprod_dim, 1)
        query_for_att = self.embed_dim_proj_query(
            query_for_att
        )  # (batch_size, outerprod_dim, n_channels * n_heads)
        kv_for_att = self.kv_proj(kv_embed).unsqueeze(-1)
        kv_for_att = self.embed_dim_proj_kv(kv_for_att)

        # Unlike with JAX, with torch we can only have one batch dimension
        # so we flatten the batch and mc samples
        if has_mc_samples:
            query_embed_flat_batch = torch.reshape(
                query_embed, (query_embed.shape[0] * query_embed.shape[1], query_embed.shape[2])
            )
            query_for_att = torch.reshape(
                query_for_att,
                (query_for_att.shape[0] * query_for_att.shape[1], query_for_att.shape[2], -1),
            )
            kv_for_att = torch.reshape(
                kv_for_att, (kv_for_att.shape[0] * kv_for_att.shape[1], kv_for_att.shape[2], -1)
            )

        eps = self.attention(query_for_att, kv_for_att, kv_for_att, need_weights=False)[
            0
        ]  # (batch_size, outerprod_dim, n_channels * n_heads)

        eps = torch.reshape(eps, (eps.shape[0], -1))
        eps_ = self.mlp_eps(eps, training)
        inputs = torch.cat(
            [query_embed_flat_batch if has_mc_samples else query_embed, eps_], dim=-1
        )
        residual = self.mlp_residual(inputs, training)

        if has_mc_samples:
            # Reshape back to the original batch and mc samples dimensions
            residual = torch.reshape(residual, (query_embed.shape[0], query_embed.shape[1], -1))

        return residual
