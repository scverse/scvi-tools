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


class _ConvBlock(nn.Module):
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
            nn.MaxPool1d(pool_size, padding=(pool_size - 1) // 2, ceil_mode=True)
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


class _DenseBlock(nn.Module):
    def __init__(
        self,
        in_features: int,
        out_features: int,
        batch_norm: bool = True,
        dropout: float = 0.2,
        activation_fn: Optional[Callable] = None,
    ):
        super().__init__()
        self.dense = nn.Linear(in_features, out_features, bias=not batch_norm)
        self.batch_norm = nn.BatchNorm1d(out_features) if batch_norm else nn.Identity()
        self.activation_fn = activation_fn if activation_fn is not None else nn.GELU()
        self.dropout = nn.Dropout(dropout)

    def forward(self, x: torch.Tensor):
        x = self.dense(x)
        x = self.batch_norm(x)
        x = self.dropout(x)
        x = self.activation_fn(x)
        return x


class _StochasticReverseComplement(nn.Module):
    """Stochastically reverse complement a one hot encoded DNA sequence."""

    def __init__(self):
        super().__init__()

    def forward(self, seq_1hot: torch.Tensor):
        """Stochastically reverse complement a one hot encoded DNA sequence.

        Parameters
        ----------
        seq_1hot
            [batch_size, seq_depth, seq_length] sequence
        """
        if self.training:
            reverse_bool = torch.rand(1).item() > 0.5
            if reverse_bool:
                # Reverse on the 4dim DNA dimension (A->T, C->G, G->C, T->A)
                # Equivalent to reversing based on our encoding
                src_seq_1hot = torch.flip(seq_1hot, [-2])
                # Reverse the sequence
                src_seq_1hot = torch.flip(src_seq_1hot, [-1])
            else:
                src_seq_1hot = seq_1hot
            return src_seq_1hot, reverse_bool
        else:
            return seq_1hot, False


class _StochasticShift(nn.Module):
    """Stochastically shift a one hot encoded DNA sequence."""

    def __init__(self, shift_max=0, pad="uniform", **kwargs):
        super().__init__()
        self.shift_max = shift_max
        self.register_buffer(
            "augment_shifts", torch.arange(-self.shift_max, self.shift_max + 1)
        )
        self.pad = pad

    def forward(self, seq_1hot: torch.tensor):
        if self.training:
            shift_i = torch.randint(0, len(self.augment_shifts), size=(1,))
            shift = self.augment_shifts[shift_i]
            if shift != 0:
                sseq_1hot = self.shift_sequence(seq_1hot, shift)
            else:
                sseq_1hot = seq_1hot
            return sseq_1hot
        else:
            return seq_1hot

    @staticmethod
    def shift_sequence(seq: torch.tensor, shift: torch.tensor, pad_value: float = 0.25):
        """Shift a sequence left or right by shift_amount.

        Parameters
        ----------
        seq
            [batch_size, seq_depth, seq_length] sequence
        shift
            signed shift value (torch.int32 or int)
        pad_value
            value to fill the padding (primitive or scalar tensor)
        """
        if len(seq.shape) != 3:
            raise ValueError("input sequence should be rank 3")
        input_shape = seq.shape

        pad = pad_value * torch.ones_like(seq[..., 0 : torch.abs(shift)])

        def _shift_right(_seq):
            # shift is positive
            sliced_seq = _seq[..., :-shift:]
            return torch.concat([pad, sliced_seq], dim=-1)

        def _shift_left(_seq):
            # shift is negative
            sliced_seq = _seq[..., -shift:]
            return torch.concat([sliced_seq, pad], dim=-1)

        sseq = _shift_right(seq) if shift > 0 else _shift_left(seq)
        sseq.reshape(input_shape)

        return sseq


class ScBassetModule(BaseModuleClass):
    """
    PyTorch implementation of ScBasset :cite:p:`Yuan2022`

    Original implementation in Keras: https://github.com/calico/scBasset

    Parameters
    ----------
    n_cells
        Number of cells to predict region accessibility
    batch_ids
        Array of (n_cells,) with batch ids for each cell
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
        batch_ids: Optional[np.ndarray] = None,
        n_filters_init: int = 288,
        n_repeat_blocks_tower: int = 6,
        filters_mult: float = 1.122,
        n_filters_pre_bottleneck: int = 256,
        n_bottleneck_layer: int = 32,
        batch_norm: bool = True,
        dropout: float = 0.0,
    ):
        super().__init__()

        self.cell_embedding = nn.Parameter(torch.randn(n_bottleneck_layer, n_cells))
        self.cell_bias = nn.Parameter(torch.randn(n_cells))
        if batch_ids is not None:
            self.register_buffer("batch_ids", torch.as_tensor(batch_ids).long())
            self.n_batch = len(torch.unique(batch_ids))
            self.register_buffer(
                "batch_ids_one_hot",
                torch.nn.functional.one_hot(batch_ids, self.n_batch)
                .squeeze(1)
                .float()
                .T,
            )
            self.batch_emdedding = nn.Parameter(
                torch.randn(n_bottleneck_layer, self.n_batch)
            )
        self.stem = _ConvBlock(
            in_channels=4,
            out_channels=n_filters_init,
            kernel_size=17,
            pool_size=3,
            dropout=dropout,
            batch_norm=batch_norm,
        )

        tower_layers = []
        curr_n_filters = n_filters_init
        for i in range(n_repeat_blocks_tower):
            new_n_filters = (
                _round(curr_n_filters * filters_mult) if i > 0 else curr_n_filters
            )
            tower_layers.append(
                _ConvBlock(
                    in_channels=curr_n_filters,
                    out_channels=new_n_filters,
                    kernel_size=5,
                    pool_size=2,
                    dropout=dropout,
                    batch_norm=batch_norm,
                )
            )
            curr_n_filters = new_n_filters
        self.tower = nn.Sequential(*tower_layers)

        self.pre_bottleneck = _ConvBlock(
            in_channels=curr_n_filters,
            out_channels=n_filters_pre_bottleneck,
            kernel_size=1,
            dropout=dropout,
            batch_norm=batch_norm,
            pool_size=1,
        )
        self.bottleneck = _DenseBlock(
            in_features=n_filters_pre_bottleneck * 7,
            out_features=n_bottleneck_layer,
            batch_norm=True,
            dropout=0.2,
            activation_fn=nn.Identity(),
        )
        self.stochastic_rc = _StochasticReverseComplement()
        self.stochastic_shift = _StochasticShift(3)

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
        h, _ = self.stochastic_rc(h)
        h = self.stochastic_shift(h)
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
        accessibility = region_embedding @ self.cell_embedding
        accessibility += self.cell_bias
        if hasattr(self, "batch_ids_one_hot"):
            accessibility += (
                region_embedding @ self.batch_emdedding
            ) @ self.batch_ids_one_hot
        return {"reconstruction": torch.sigmoid(accessibility)}

    def loss(self, tensors, inference_outputs, generative_outputs) -> LossOutput:
        """Loss function for the model."""
        reconstruction = generative_outputs["reconstruction"]
        target = tensors[REGISTRY_KEYS.X_KEY]
        loss_fn = nn.BCELoss()
        loss = loss_fn(reconstruction, target)
        return LossOutput(loss=loss, n_obs_minibatch=target.shape[0])
