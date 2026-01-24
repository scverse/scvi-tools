"""Base neural network components for DIAGVI model."""

from __future__ import annotations

import torch
import torch.nn as nn
import torch.nn.functional as F
from torch.distributions import Normal

from scvi.utils import dependencies

EPS = 1e-8


class DecoderRNA(nn.Module):
    """Decoder for RNA modality using feature embeddings.

    Decodes latent representations to RNA expression using batch-specific
    scale and bias parameters combined with feature embeddings from the
    guidance graph.

    Parameters
    ----------
    n_output
        Number of output features (genes).
    n_batches
        Number of batches in the data.
    """

    def __init__(
        self,
        n_output: int,
        n_batches: int,
    ):
        super().__init__()
        self.n_output = n_output
        self.n_batches = n_batches

        self.scale_lin = nn.Parameter(torch.zeros(n_batches, n_output))
        self.bias = nn.Parameter(torch.zeros(n_batches, n_output))
        self.log_theta = nn.Parameter(torch.zeros(n_batches, n_output))

        self.px_dropout_param = nn.Parameter(torch.randn(n_output) * 0.01)

    def forward(
        self,
        u: torch.Tensor,
        l: torch.Tensor,
        batch_index: torch.Tensor,
        v: torch.Tensor,
    ) -> tuple[torch.Tensor, torch.Tensor, torch.Tensor, torch.Tensor]:
        """Decode latent representation to RNA expression parameters.

        Parameters
        ----------
        u
            Latent representation tensor of shape (n_cells, n_latent).
        l
            Log library size tensor of shape (n_cells, 1).
        batch_index
            Batch indices tensor of shape (n_cells,) or (n_cells, 1).
        v
            Feature embedding tensor from graph encoder.

        Returns
        -------
        tuple[torch.Tensor, torch.Tensor, torch.Tensor, torch.Tensor]
            Tuple of (px_scale, px_r, px_rate, px_dropout) for negative binomial
            distribution parameters.
        """
        if batch_index.dim() > 1:
            batch_index = batch_index.squeeze(-1)

        if (batch_index.max() >= self.bias.shape[0]) or (batch_index.min() < 0):
            raise IndexError(
                f"Batch index out of bounds: valid range is [0, {self.bias.shape[0] - 1}]"
            )

        scale = F.softplus(self.scale_lin[batch_index])
        bias = self.bias[batch_index]
        log_theta = self.log_theta[batch_index]

        raw_px_scale = scale * (u @ v.T) + bias
        px_scale = torch.softmax(raw_px_scale, dim=-1)
        px_rate = torch.exp(l) * px_scale

        px_dropout = F.softplus(self.px_dropout_param)

        px_r = log_theta

        return px_scale, px_r, px_rate, px_dropout


class DecoderProteinGLUE(nn.Module):
    """Decoder for protein modality using GLUE-style mixture model.

    Decodes latent representations to protein expression using a mixture
    of two components with batch-specific parameters.

    Parameters
    ----------
    n_output
        Number of output features (proteins).
    n_batches
        Number of batches in the data.
    """

    def __init__(
        self,
        n_output: int,
        n_batches: int,
    ):
        super().__init__()
        self.n_output = n_output
        self.n_batches = n_batches

        self.scale_lin = nn.Parameter(torch.zeros(n_batches, n_output))

        self.bias1 = nn.Parameter(torch.zeros(n_batches, n_output))
        self.bias2 = nn.Parameter(torch.zeros(n_batches, n_output))

        self.log_theta = nn.Parameter(torch.zeros(n_batches, n_output))

    def forward(
        self,
        u: torch.Tensor,
        l: torch.Tensor,
        batch_index: torch.Tensor,
        v: torch.Tensor,
    ) -> tuple[
        tuple[torch.Tensor, torch.Tensor],
        torch.Tensor,
        tuple[torch.Tensor, torch.Tensor],
        torch.Tensor,
    ]:
        """Decode latent representation to protein expression parameters.

        Parameters
        ----------
        u
            Latent representation tensor of shape (n_cells, n_latent).
        l
            Log library size tensor of shape (n_cells, 1).
        batch_index
            Batch indices tensor of shape (n_cells,) or (n_cells, 1).
        v
            Feature embedding tensor from graph encoder.

        Returns
        -------
        tuple
            Tuple of ((px_scale_1, px_scale_2), px_r, (px_rate_1, px_rate_2),
            mixture_logits) for negative binomial mixture distribution.
        """
        if batch_index.dim() > 1:
            batch_index = batch_index.squeeze(-1)

        if (batch_index.max() >= self.bias1.shape[0]) or (batch_index.min() < 0):
            raise IndexError(
                f"Batch index out of bounds: valid range is [0, {self.bias1.shape[0] - 1}]"
            )

        scale = F.softplus(self.scale_lin[batch_index])

        bias1 = self.bias1[batch_index]
        bias2 = self.bias2[batch_index]

        log_theta = self.log_theta[batch_index]

        raw_px_scale_1 = scale * (u @ v.T) + bias1
        raw_px_scale_2 = scale * (u @ v.T) + bias2

        px_scale_1 = torch.softmax(raw_px_scale_1, dim=-1)
        px_scale_2 = torch.softmax(raw_px_scale_2, dim=-1)

        px_rate_1 = torch.exp(l) * px_scale_1
        px_rate_2 = torch.exp(l) * px_scale_2

        mixture_logits = raw_px_scale_1 - raw_px_scale_2

        px_r = log_theta

        return (px_scale_1, px_scale_2), px_r, (px_rate_1, px_rate_2), mixture_logits


class GraphEncoder(nn.Module):
    """Graph convolutional encoder for feature embeddings.

    Encodes feature nodes in the guidance graph to learn feature
    embeddings that capture cross-modality relationships.

    Parameters
    ----------
    vnum
        Number of nodes (features) in the graph.
    out_features
        Dimensionality of the output feature embeddings.
    """

    @dependencies("torch_geometric")
    def __init__(self, vnum: int, out_features: int):
        import torch_geometric

        super().__init__()
        self.vrepr = nn.Parameter(torch.zeros(vnum, out_features))
        self.conv = torch_geometric.nn.GCNConv(out_features, out_features)
        self.loc = nn.Linear(out_features, out_features)
        self.std_lin = nn.Linear(out_features, out_features)

    def forward(self, edge_index: torch.Tensor) -> tuple[torch.Tensor, torch.Tensor, torch.Tensor]:
        """Encode graph to feature embeddings.

        Parameters
        ----------
        edge_index
            Edge index tensor of shape (2, n_edges).

        Returns
        -------
        tuple[torch.Tensor, torch.Tensor, torch.Tensor]
            Tuple of (z, mu, logvar) where z is the sampled embedding,
            mu is the mean, and logvar is the log variance.
        """
        h = self.conv(self.vrepr, edge_index)
        loc = self.loc(h)
        std = F.softplus(self.std_lin(h)) + EPS

        dist = Normal(loc, std)
        mu = dist.loc
        std = dist.scale
        logvar = torch.log(std**2)
        z = dist.rsample()

        return z, mu, logvar
