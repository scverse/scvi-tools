import torch
import torch.nn as nn
import torch.nn.functional as F
from torch.distributions import Normal
from torch_geometric.nn import GCNConv

EPS = 1e-8


class NBDataDecoderWB(nn.Module):  # integrate the batch index
    def __init__(
        self,
        n_output: int,
        n_latent: int,
        n_batches: int,
    ):
        super().__init__()
        self.n_output = n_output
        self.n_batches = n_batches

        self.scale_lin = nn.Parameter(torch.zeros(n_batches, n_output))
        self.bias = nn.Parameter(torch.zeros(n_batches, n_output))
        self.log_theta = nn.Parameter(torch.zeros(n_batches, n_output))

        self.px_dropout_param = nn.Parameter(torch.randn(n_output) * 0.01)
        # self.v = nn.Parameter(torch.randn(n_output, n_latent) * 0.01)

    def forward(
        self,
        u: torch.Tensor,
        l: torch.Tensor,
        batch_index: torch.Tensor,
        v: torch.Tensor,
    ) -> tuple[torch.Tensor, torch.Tensor]:
        # bring batch index in the right dimension
        if batch_index.dim() > 1:
            batch_index = batch_index.squeeze(-1)

        scale = F.softplus(self.scale_lin[batch_index])

        if (batch_index.max() >= self.bias.shape[0]) or (batch_index.min() < 0):
            raise IndexError(
                f"Batch index out of bounds: valid range is [0, {self.bias.shape[0] - 1}]"
            )

        bias = self.bias[batch_index]
        log_theta = self.log_theta[batch_index]

        raw_px_scale = scale * (u @ v.T) + bias
        px_scale = torch.softmax(raw_px_scale, dim=-1)
        px_rate = torch.exp(l) * px_scale

        px_dropout = F.softplus(self.px_dropout_param)

        px_r = log_theta

        return px_scale, px_r, px_rate, px_dropout


class GraphEncoder_glue(nn.Module):
    def __init__(self, vnum: int, out_features: int):
        super().__init__()
        self.vrepr = nn.Parameter(torch.zeros(vnum, out_features))
        self.conv = GCNConv(out_features, out_features)  # evtl auch GAT - user parameter
        self.loc = nn.Linear(out_features, out_features)
        self.std_lin = nn.Linear(out_features, out_features)

    def forward(self, edge_index):
        h = self.conv(self.vrepr, edge_index)
        loc = self.loc(h)
        std = F.softplus(self.std_lin(h)) + EPS

        dist = Normal(loc, std)
        mu = dist.loc
        std = dist.scale
        logvar = torch.log(std**2)
        z = dist.rsample()

        return z, mu, logvar
