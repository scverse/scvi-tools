from typing import Callable, Dict, NamedTuple, Optional

import numpy as np
import torch
from torch import nn

from scvi.module.base import BaseModuleClass, LossOutput, auto_move_data


class _REGISTRY_KEYS_NT(NamedTuple):
    X_KEY: str = "X"
    BATCH_KEY: str = "batch"
    DNA_CODE_KEY: str = "dna_code"


REGISTRY_KEYS = _REGISTRY_KEYS_NT()


def _round(x):
    return int(np.round(x))


class _ConvLayer(nn.Module):
    def __init__(
        self,
        in_channels: int,
        out_channels: int,
        kernel_size: int,
        pool_size: int = None,
        batch_norm: bool = True,
        dropout: float = 0.0,
        activation_fn: Optional[Callable] = None,
    ):
        super().__init__()
        self.conv = nn.Conv1d(
            in_channels=in_channels,
            out_channels=out_channels,
            kernel_size=kernel_size,
            padding="same",
            bias=False,
        )
        self.batch_norm = nn.BatchNorm1d(out_channels) if batch_norm else nn.Identity()
        self.pool = (
            nn.MaxPool1d(pool_size, padding=(pool_size - 1) // 2)
            if pool_size is not None
            else nn.Identity()
        )
        self.activation_fn = activation_fn if activation_fn is not None else nn.GELU()
        self.dropout = nn.Dropout(dropout)

    def forward(self, x: torch.Tensor):
        x = self.conv(x)
        x = self.batch_norm(x)
        x = self.dropout(x)
        x = self.pool(x)
        x = self.activation_fn(x)
        return x


class _DenseLayer(nn.Module):
    def __init__(
        self,
        in_features: int,
        out_features: int,
        use_bias: bool = True,
        batch_norm: bool = True,
        dropout: float = 0.2,
        activation_fn: Optional[Callable] = None,
    ):
        super().__init__()
        self.dense = nn.Linear(in_features, out_features, bias=use_bias)
        self.batch_norm = nn.BatchNorm1d(out_features) if batch_norm else nn.Identity()
        self.activation_fn = activation_fn if activation_fn is not None else nn.GELU()
        self.dropout = nn.Dropout(dropout)

    def forward(self, x: torch.Tensor):
        x = self.dense(x)
        x = self.batch_norm(x)
        x = self.dropout(x)
        x = self.activation_fn(x)
        return x


class ScBassetModule(BaseModuleClass):
    """
    PyTorch implementation of ScBasset :cite:p:`Yuan2022`

    Original implementation in Keras: https://github.com/calico/scBasset

    Parameters
    ----------
    n_cells
        Number of cells to predict region accessibility
    n_filters_init
        Number of filters for the initial conv layer
    n_repeat_blocks_tower
        Number of layers in the convolutional tower
    filters_mult
        Proportion by which the number of filters should inrease in the
        convolutional tower
    n_bottleneck_layer
        Size of the bottleneck layer
    batch_norm
        Whether to apply batch norm across model layers
    dropout
        Dropout rate across layers, by default we do not do it for
        convolutional layers but we do it for the dense layers
    """

    def __init__(
        self,
        n_cells: int,
        n_filters_init: int = 288,
        n_repeat_blocks_tower: int = 5,
        filters_mult: float = 1.122,
        n_filters_pre_bottleneck: int = 256,
        n_bottleneck_layer: int = 32,
        batch_norm: bool = True,
        dropout: float = 0.0,
    ):
        super().__init__()

        self.cell_embedding = nn.Parameter(torch.randn(n_cells, n_bottleneck_layer))
        self.stem = _ConvLayer(
            in_channels=4,
            out_channels=n_filters_init,
            kernel_size=17,
            pool_size=3,
            dropout=dropout,
            batch_norm=batch_norm,
        )

        tower_layers = []
        curr_n_filters = n_filters_init
        for _ in range(n_repeat_blocks_tower):
            tower_layers.append(
                _ConvLayer(
                    in_channels=curr_n_filters,
                    out_channels=_round(curr_n_filters * filters_mult),
                    kernel_size=5,
                    pool_size=2,
                    dropout=dropout,
                    batch_norm=batch_norm,
                )
            )
            curr_n_filters = _round(curr_n_filters * filters_mult)
        self.tower = nn.Sequential(*tower_layers)

        self.pre_bottleneck = _ConvLayer(
            in_channels=curr_n_filters,
            out_channels=n_filters_pre_bottleneck,
            kernel_size=1,
            dropout=dropout,
            batch_norm=batch_norm,
            pool_size=2,
        )
        self.bottleneck = _DenseLayer(
            in_features=n_filters_pre_bottleneck * 6,
            out_features=n_bottleneck_layer,
            use_bias=True,
            batch_norm=True,
            dropout=0.2,
            activation_fn=nn.Identity(),
        )

    def _get_inference_input(self, tensors: Dict[str, torch.Tensor]):
        dna_code = tensors[REGISTRY_KEYS.DNA_CODE_KEY]

        input_dict = dict(dna_code=dna_code)
        return input_dict

    @auto_move_data
    def inference(self, dna_code: torch.Tensor) -> Dict[str, torch.Tensor]:
        """Inference method for the model."""
        # input shape: (batch_size, seq_length)
        # output shape: (batch_size, 4, seq_length)
        h = nn.functional.one_hot(dna_code, num_classes=4).permute(0, 2, 1).float()

        # TODO: add random shift to act as a regularizer on the dataset level
        # TODO: add use reverse complement randomly on the dataset level
        h = self.stem(h)
        h = self.tower(h)
        h = self.pre_bottleneck(h)
        # flatten the input
        h = h.view(h.shape[0], -1)
        # Regions by bottleneck layer dim
        h = self.bottleneck(h)
        return {"region_embedding": h}

    def _get_generative_input(
        self,
        tensors: Dict[str, torch.Tensor],
        inference_outputs: Dict[str, torch.Tensor],
    ):
        region_embedding = inference_outputs["region_embedding"]
        input_dict = dict(region_embedding=region_embedding)

        return input_dict

    def generative(self, region_embedding: torch.Tensor) -> Dict[str, torch.Tensor]:
        """Generative method for the model."""
        return {
            "reconstruction": nn.functional.sigmoid(
                region_embedding @ self.cell_embedding.T
            )
        }

    def loss(self, tensors, inference_outputs, generative_outputs) -> LossOutput:
        """Loss function for the model."""
        reconstruction = generative_outputs["reconstruction"]
        target = tensors[REGISTRY_KEYS.X_KEY]
        loss = nn.functional.binary_cross_entropy(reconstruction, target)
        return LossOutput(loss=loss, n_obs_minibatch=target.shape[0])
