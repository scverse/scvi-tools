"""Base neural network components for DIAGVI model."""

from __future__ import annotations

import torch
import torch.nn as nn
import torch.nn.functional as F
from torch.distributions import Normal

from scvi.utils import dependencies

EPS = 1e-8


class DecoderSinglePathway(nn.Module):
    """Single-pathway decoder for unimodal likelihoods.

    Outputs a single set of parameters (scale, r, rate, dropout).
    Use for: nb, zinb, normal, log1pnormal, ziln, zig.

    Parameters
    ----------
    n_output
        Number of output features.
    n_batches
        Number of batches.
    normalize
        If True, apply softmax normalization to output scale (for count data).
        If False, output raw scale values (for continuous data).
        Default is True for backward compatibility.
    """

    def __init__(
        self,
        n_output: int,
        n_batches: int,
        normalize: bool = True,
    ):
        super().__init__()
        self.n_output = n_output
        self.n_batches = n_batches
        self.normalize = normalize

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
        """Forward pass through the decoder.

        Parameters
        ----------
        u
            Cell latent representation, shape (n_cells, n_latent).
        l
            Library size (log-scale), shape (n_cells, 1).
        batch_index
            Batch indices for cells, shape (n_cells,) or (n_cells, 1).
        v
            Feature latent representation, shape (n_features, n_latent).

        Returns
        -------
        tuple of 4 :class:`torch.Tensor`
            - px_scale: Normalized scale parameters, shape (n_cells, n_features).
            - px_r: Overdispersion parameter, shape (n_cells, n_features).
            - px_rate: Rate parameters (mean), shape (n_cells, n_features).
            - px_dropout: Dropout/zero-inflation probabilities, shape (n_features,).
        """
        # Check batch index dimension
        if batch_index.dim() > 1:
            batch_index = batch_index.squeeze(-1)

        if (batch_index.max() >= self.bias.shape[0]) or (batch_index.min() < 0):
            raise IndexError(
                f"Batch index out of bounds: valid range is [0, {self.bias.shape[0] - 1}]"
            )

        # Get batch-specific parameters
        scale = F.softplus(self.scale_lin[batch_index])
        bias = self.bias[batch_index]
        log_theta = self.log_theta[batch_index]

        # Compute decoder output parameters
        raw_px_scale = scale * (u @ v.T) + bias

        if self.normalize:
            # Count data: softmax produces proportions, library scales to counts
            px_scale = torch.softmax(raw_px_scale, dim=-1)
            px_rate = torch.exp(l) * px_scale
        else:
            # Continuous data: no normalization, scale is positive, rate is unbounded
            px_scale = F.softplus(raw_px_scale) + EPS
            px_rate = raw_px_scale

        px_dropout = F.softplus(self.px_dropout_param)
        px_r = log_theta

        return px_scale, px_r, px_rate, px_dropout


class DecoderDualPathway(nn.Module):
    """Dual-pathway decoder for mixture likelihoods.

    Outputs two sets of parameters for foreground/background mixture models.
    Use for: nbmixture.

    Parameters
    ----------
    n_output
        Number of output features.
    n_batches
        Number of batches.
    normalize
        If True, apply softmax normalization to output scales (for count data).
        If False, output raw scale values (for continuous data).
        Default is True for backward compatibility.
    """

    def __init__(
        self,
        n_output: int,
        n_batches: int,
        normalize: bool = True,
    ):
        super().__init__()
        self.n_output = n_output
        self.n_batches = n_batches
        self.normalize = normalize

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
    ) -> tuple[tuple[torch.Tensor, torch.Tensor], torch.Tensor, tuple[torch.Tensor, torch.Tensor], torch.Tensor]:
        """Forward pass through the decoder.

        Parameters
        ----------
        u
            Cell latent representation, shape (n_cells, n_latent).
        l
            Library size (log-scale), shape (n_cells, 1).
        batch_index
            Batch indices for cells, shape (n_cells,) or (n_cells, 1).
        v
            Feature latent representation, shape (n_features, n_latent).

        Returns
        -------
        tuple of 4 elements
            - scales: Tuple of (foreground, background) scale parameters.
              Each has shape (n_cells, n_features).
            - px_r: Overdispersion parameter, shape (n_cells, n_features).
            - rates: Tuple of (foreground, background) rate parameters (mean).
              Each has shape (n_cells, n_features).
            - mixture_logits: Logits for mixture component probabilities, shape (n_cells, n_features).
        """
        # Check batch index dimension
        if batch_index.dim() > 1:
            batch_index = batch_index.squeeze(-1)

        if (batch_index.max() >= self.bias1.shape[0]) or (batch_index.min() < 0):
            raise IndexError(
                f"Batch index out of bounds: valid range is [0, {self.bias1.shape[0] - 1}]"
            )

        # Get batch-specific parameters
        scale = F.softplus(self.scale_lin[batch_index])
        bias1 = self.bias1[batch_index]
        bias2 = self.bias2[batch_index]
        log_theta = self.log_theta[batch_index]

        # Compute decoder output parameters for both pathways
        raw_px_scale_1 = scale * (u @ v.T) + bias1
        raw_px_scale_2 = scale * (u @ v.T) + bias2

        if self.normalize:
            # Count data: softmax produces proportions, library scales to counts
            px_scale_1 = torch.softmax(raw_px_scale_1, dim=-1)
            px_scale_2 = torch.softmax(raw_px_scale_2, dim=-1)
            px_rate_1 = torch.exp(l) * px_scale_1
            px_rate_2 = torch.exp(l) * px_scale_2
        else:
            # Continuous data: no normalization, scales are positive, rates are unbounded
            px_scale_1 = F.softplus(raw_px_scale_1) + EPS
            px_scale_2 = F.softplus(raw_px_scale_2) + EPS
            px_rate_1 = raw_px_scale_1
            px_rate_2 = raw_px_scale_2

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
        """Forward pass through the graph encoder.

        Parameters
        ----------
        edge_index
            Edge indices for the graph, shape (2, n_edges).

        Returns
        -------
        tuple of 3 :class:`torch.Tensor`
            - z: Sampled latent representation from the normal distribution, shape (n_nodes, out_features).
            - mu: Mean of the latent distribution, shape (n_nodes, out_features).
            - logvar: Log variance of the latent distribution, shape (n_nodes, out_features).
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
