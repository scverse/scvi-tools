import logging
from typing import Callable, NamedTuple, Optional

import numpy as np
import torch
import torchmetrics
from torch import nn

from scvi.module.base import BaseModuleClass, LossOutput, auto_move_data

logger = logging.getLogger(__name__)


class _REGISTRY_KEYS_NT(NamedTuple):
    X_KEY: str = "X"
    BATCH_KEY: str = "batch"
    DNA_CODE_KEY: str = "dna_code"


REGISTRY_KEYS = _REGISTRY_KEYS_NT()


def _round(x):
    return int(np.round(x))


class _Linear(nn.Linear):
    """Linear layer with Keras default initalizations."""

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

    def reset_parameters(self):
        """Reset parameters."""
        torch.nn.init.kaiming_normal_(self.weight)
        if self.bias is not None:
            torch.nn.init.zeros_(self.bias)


class _GELU(nn.Module):
    """GELU unit approximated by a sigmoid, same as original."""

    def __init__(self):
        super().__init__()

    def forward(self, x: torch.Tensor):
        return torch.sigmoid(1.702 * x) * x


class _BatchNorm(nn.BatchNorm1d):
    """Batch normalization with Keras default initializations and scBasset defaults."""

    def __init__(self, *args, eps: int = 1e-3, **kwargs):
        # Keras uses 0.01 for momentum, but scBasset uses 0.1
        super().__init__(*args, eps=eps, **kwargs)


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
        ceil_mode: bool = False,
    ):
        super().__init__()
        self.conv = nn.Conv1d(
            in_channels=in_channels,
            out_channels=out_channels,
            kernel_size=kernel_size,
            padding="same",
            bias=False,
        )
        self.batch_norm = _BatchNorm(out_channels) if batch_norm else nn.Identity()
        self.pool = (
            nn.MaxPool1d(pool_size, padding=(pool_size - 1) // 2, ceil_mode=ceil_mode)
            if pool_size is not None
            else nn.Identity()
        )
        self.activation_fn = activation_fn if activation_fn is not None else _GELU()
        self.dropout = nn.Dropout(dropout)

    def forward(self, x: torch.Tensor):
        x = self.activation_fn(x)
        x = self.conv(x)
        x = self.batch_norm(x)
        x = self.dropout(x)
        x = self.pool(x)
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
        self.dense = _Linear(in_features, out_features, bias=not batch_norm)
        # batch norm with Keras default epsilon
        self.batch_norm = _BatchNorm(out_features) if batch_norm else nn.Identity()
        self.activation_fn = activation_fn if activation_fn is not None else _GELU()
        self.dropout = nn.Dropout(dropout)

    def forward(self, x: torch.Tensor):
        x = self.activation_fn(x)
        x = self.dense(x)
        x = self.batch_norm(x)
        x = self.dropout(x)
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
            reverse_bool = np.random.uniform() > 0.5
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
        self.augment_shifts = np.arange(-self.shift_max, self.shift_max + 1)
        self.pad = pad

    def forward(self, seq_1hot: torch.Tensor):
        if self.training:
            shift_i = np.random.randint(0, len(self.augment_shifts))
            shift = self.augment_shifts[shift_i]
            if shift != 0:
                return self.shift_sequence(seq_1hot, shift)
            else:
                return seq_1hot
        else:
            return seq_1hot

    @staticmethod
    def shift_sequence(seq: torch.Tensor, shift: int, pad_value: float = 0.25):
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

        sseq = torch.roll(seq, shift, dims=-1)
        if shift > 0:
            sseq[..., :shift] = pad_value
        else:
            sseq[..., shift:] = pad_value

        return sseq


class ScBassetModule(BaseModuleClass):
    """PyTorch implementation of ScBasset :cite:p:`Yuan2022`.

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
    l2_reg_cell_embedding
        L2 regularization for the cell embedding layer
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
        l2_reg_cell_embedding: float = 0.0,
    ):
        super().__init__()
        self.l2_reg_cell_embedding = l2_reg_cell_embedding
        self.cell_embedding = nn.Parameter(torch.randn(n_bottleneck_layer, n_cells))
        self.cell_bias = nn.Parameter(torch.randn(n_cells))
        if batch_ids is not None:
            self.register_buffer("batch_ids", torch.as_tensor(batch_ids).long())
            self.n_batch = len(torch.unique(batch_ids))
            self.batch_emdedding = nn.Embedding(self.n_batch, n_bottleneck_layer)
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
        # NOTE: Bottleneck here assumes that seq_len=1344 and n_repeat_blocks_tower=6
        # seq_len and tower size are fixed by the in_features shape
        self.bottleneck = _DenseBlock(
            in_features=n_filters_pre_bottleneck * 7,
            out_features=n_bottleneck_layer,
            batch_norm=True,
            dropout=0.2,
        )
        self.stochastic_rc = _StochasticReverseComplement()
        self.stochastic_shift = _StochasticShift(3)

    def _get_inference_input(self, tensors: dict[str, torch.Tensor]):
        dna_code = tensors[REGISTRY_KEYS.DNA_CODE_KEY]

        input_dict = {"dna_code": dna_code}
        return input_dict

    @auto_move_data
    def inference(self, dna_code: torch.Tensor) -> dict[str, torch.Tensor]:
        """Inference method for the model."""
        # NOTE: `seq_len` assumed to be a fixed 1344 as in the original implementation.
        # input shape: (batch_size, seq_length)
        # output shape: (batch_size, 4, seq_length)
        h = nn.functional.one_hot(dna_code, num_classes=4).permute(0, 2, 1).float()
        h, _ = self.stochastic_rc(h)
        h = self.stochastic_shift(h)
        # input shape: (batch_size, 4, seq_length)
        # output shape: (batch_size, n_filters_stem, seq_length//3)
        # `stem` contains a max_pool1d by 3. For 1344 input, now 448
        h = self.stem(h)
        # output shape: (batch_size, n_filters_tower, seq_length//(3*2**6))
        # `tower` contains 6 (default) `max_pool1d` by 2
        # for 1344 input, now 7
        h = self.tower(h)
        # output shape: (batch_size, n_filters_pre_bottleneck=1, seq_length//(3*2**6))
        # `bottleneck` is a filter k=1 conv with no pooling
        # for 1344 input, now [batch_size, 1, 7]
        h = self.pre_bottleneck(h)
        # flatten the input
        # output shape: (batch_size, n_filters_pre_bottleneck * (seq_length//(3*2**6)))
        # for 1344 input, now [batch_size, 7]
        h = h.view(h.shape[0], -1)
        # Regions by bottleneck layer dim
        # output shape: (batch_size, n_bottleneck_layer)
        h = self.bottleneck(h)
        h = _GELU()(h)
        return {"region_embedding": h}

    def _get_generative_input(
        self,
        tensors: dict[str, torch.Tensor],
        inference_outputs: dict[str, torch.Tensor],
    ):
        region_embedding = inference_outputs["region_embedding"]
        input_dict = {"region_embedding": region_embedding}

        return input_dict

    def _get_accessibility(
        self,
        dna_codes: torch.Tensor,
        batch_size: int = None,
    ) -> torch.Tensor:
        """Perform minibatch inference of accessibility scores."""
        accessibility = torch.zeros(
            size=(
                dna_codes.shape[0],
                self.cell_bias.shape[0],
            )
        )
        if batch_size is None:
            # no minibatching
            batch_size = dna_codes.shape[0]
        n_batches = accessibility.shape[0] // batch_size + int(
            (accessibility.shape[0] % batch_size) > 0
        )
        for batch in range(n_batches):
            batch_codes = dna_codes[batch * batch_size : (batch + 1) * batch_size]
            # forward passes, output is dict(region_embedding=np.ndarray: [n_seqs, n_latent=32])
            motif_rep = self.inference(dna_code=batch_codes)
            # output is dict(reconstruction_logits=np.ndarray: [n_seqs, n_cells])
            batch_acc = self.generative(region_embedding=motif_rep["region_embedding"])[
                "reconstruction_logits"
            ]
            accessibility[batch * batch_size : (batch + 1) * batch_size] = batch_acc
        return accessibility

    def generative(self, region_embedding: torch.Tensor) -> dict[str, torch.Tensor]:
        """Generative method for the model."""
        if hasattr(self, "batch_ids"):
            # embeddings dim by cells dim
            cell_batch_embedding = self.batch_emdedding(self.batch_ids).squeeze(-2).T
        else:
            cell_batch_embedding = 0
        accessibility = region_embedding @ (self.cell_embedding + cell_batch_embedding)
        accessibility += self.cell_bias
        return {"reconstruction_logits": accessibility}

    def loss(self, tensors, inference_outputs, generative_outputs) -> LossOutput:
        """Loss function for the model."""
        reconstruction_logits = generative_outputs["reconstruction_logits"]
        target = tensors[REGISTRY_KEYS.X_KEY]
        loss_fn = nn.BCEWithLogitsLoss(reduction="none")
        full_loss = loss_fn(reconstruction_logits, target)
        reconstruction_loss = full_loss.sum(dim=-1)
        loss = reconstruction_loss.sum() / (
            reconstruction_logits.shape[0] * reconstruction_logits.shape[1]
        )
        if self.l2_reg_cell_embedding > 0:
            loss += self.l2_reg_cell_embedding * torch.square(self.cell_embedding).mean()
        auroc = torchmetrics.functional.auroc(
            torch.sigmoid(reconstruction_logits).ravel(),
            target.int().ravel(),
            task="binary",
        )
        return LossOutput(
            loss=loss,
            reconstruction_loss=reconstruction_loss,
            extra_metrics={"auroc": auroc},
        )
