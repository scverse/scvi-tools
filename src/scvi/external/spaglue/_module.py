import torch
import torch.nn as nn
import torch.nn.functional as F
from torch.distributions import Normal, kl_divergence
from torch_geometric.nn import GCNConv

from scvi import REGISTRY_KEYS
from scvi.distributions import NegativeBinomial, ZeroInflatedNegativeBinomial
from scvi.module._constants import MODULE_KEYS
from scvi.module.base import BaseModuleClass, LossOutput, auto_move_data
from scvi.nn import Encoder

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
        bias = self.bias[batch_index]
        log_theta = self.log_theta[batch_index]

        raw_px_scale = scale * (u @ v.T) + bias
        px_scale = torch.softmax(raw_px_scale, dim=-1)
        px_rate = torch.exp(l) * px_scale

        px_dropout = F.softplus(self.px_dropout_param)

        px_r = log_theta

        return px_scale, px_r, px_rate, px_dropout


class GraphEncoder(nn.Module):
    def __init__(self, in_dim, hidden_dim, out_dim):
        super().__init__()

        # GCN layers for computing mu and logvar
        self.gcn_mu_1 = GCNConv(in_dim, hidden_dim)
        self.gcn_mu_2 = GCNConv(hidden_dim, out_dim)

        self.gcn_logvar_1 = GCNConv(in_dim, hidden_dim)
        self.gcn_logvar_2 = GCNConv(hidden_dim, out_dim)

    def forward(self, x, edge_index):
        # First pass for mean
        mu = torch.relu(self.gcn_mu_1(x, edge_index))
        mu = self.gcn_mu_2(mu, edge_index)

        # Second pass for log-variance
        logvar = torch.relu(self.gcn_logvar_1(x, edge_index))
        logvar = self.gcn_logvar_2(logvar, edge_index)

        # Reparameterization trick
        std = torch.exp(0.5 * logvar)
        eps = torch.randn_like(std)
        z = mu + eps * std

        return z, mu, logvar


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


class SPAGLUEVAE(BaseModuleClass):
    def __init__(
        self,
        n_inputs: list[int],
        n_batches: list[int],
        gene_likelihoods: list[str],
        guidance_graph,
        n_latent_seq: int = 50,
        n_latent_spatial: int = 50,
        n_hidden: int = 256,
        n_layers: int = 2,
        dropout_rate: float = 0.1,
        # **kwargs: dict,
    ) -> None:
        super().__init__()

        self.n_input_list = n_inputs
        self.n_batches_list = n_batches
        self.gene_likelihoods = gene_likelihoods
        self.guidance_graph = guidance_graph

        self.z_encoder_diss = Encoder(
            n_input=n_inputs[0],
            n_output=n_latent_seq,
            n_hidden=n_hidden,
            n_layers=n_layers,
            dropout_rate=dropout_rate,
            return_dist=True,
            # **kwargs,
        )

        self.z_encoder_spa = Encoder(
            n_input=n_inputs[1],
            n_output=n_latent_spatial,
            n_hidden=n_hidden,
            n_layers=n_layers,
            dropout_rate=dropout_rate,
            return_dist=True,
            # **kwargs,
        )

        self.z_decoder_diss = NBDataDecoderWB(
            n_output=n_inputs[0],
            n_latent=n_latent_seq,
            n_batches=n_batches[0],
        )

        self.z_decoder_spa = NBDataDecoderWB(
            n_output=n_inputs[1],
            n_latent=n_latent_spatial,
            n_batches=n_batches[1],
        )

        self.graph_encoder = GraphEncoder_glue(
            vnum=n_inputs[0] + n_inputs[1],
            out_features=50,
        )

        # self.graph_encoder = GraphEncoder(
        #     in_dim=n_inputs[0] + n_inputs[1],
        #     hidden_dim=32,
        #     out_dim=10,
        # )

    def _get_inference_input(
        self, tensors: dict[str, torch.Tensor]
    ) -> dict[str, torch.Tensor | None]:
        return {
            MODULE_KEYS.X_KEY: tensors[REGISTRY_KEYS.X_KEY],
        }

    def _get_generative_input(
        self, tensors: dict[str, torch.Tensor], inference_outputs: dict[str, torch.Tensor]
    ) -> dict[str, torch.Tensor]:
        return {
            MODULE_KEYS.Z_KEY: inference_outputs[MODULE_KEYS.Z_KEY],
            MODULE_KEYS.LIBRARY_KEY: inference_outputs[MODULE_KEYS.LIBRARY_KEY],
            MODULE_KEYS.BATCH_INDEX_KEY: tensors[REGISTRY_KEYS.BATCH_KEY],
            MODULE_KEYS.Y_KEY: tensors[REGISTRY_KEYS.LABELS_KEY],
            "v": inference_outputs["v"],
        }

    @auto_move_data
    def inference(
        self,
        x: torch.Tensor,
        mode: int | None = 0,
    ) -> dict[str, torch.Tensor]:
        x_ = x
        library = torch.log(x.sum(1)).unsqueeze(1)

        graph = self.guidance_graph
        device = x.device
        graph = graph.to(device)

        # whole embedding is calculated
        v_all, mu_all, logvar_all = self.graph_encoder(graph.edge_index)

        # embedding for modality is extracted to be used for decoder input
        if mode == 0:
            v = v_all[graph.seq_indices]
        elif mode == 1:
            v = v_all[graph.spa_indices]
        else:
            raise ValueError("Invalid mode: must be 0 or 1.")

        # diss data
        if mode == 0:
            qz, z = self.z_encoder_diss(x_)
        # spa data
        else:
            qz, z = self.z_encoder_spa(x_)

        return {
            MODULE_KEYS.QZ_KEY: qz,
            MODULE_KEYS.Z_KEY: z,
            MODULE_KEYS.LIBRARY_KEY: library,
            "v": v,
            "v_all": v_all,
            "mu_all": mu_all,
            "logvar_all": logvar_all,
        }

    @auto_move_data
    def generative(
        self,
        z: torch.Tensor,
        library: torch.Tensor,
        batch_index: torch.Tensor | None = None,
        y: torch.Tensor | None = None,
        v: torch.Tensor | None = None,
        mode: int | None = 0,
    ) -> dict[str, torch.Tensor]:
        """Run the generative model."""
        EPS = 1e-8

        # diss data
        if mode == 0:
            px_scale, px_r, px_rate, px_dropout = self.z_decoder_diss(z, library, batch_index, v)

        # spa data
        else:
            px_scale, px_r, px_rate, px_dropout = self.z_decoder_spa(z, library, batch_index, v)

        px_r = px_r.exp()

        if self.gene_likelihoods[mode] == "nb":
            px = NegativeBinomial(px_r, logits=(px_rate + EPS).log() - px_r)

        elif self.gene_likelihoods[mode] == "zinb":
            px = ZeroInflatedNegativeBinomial(
                mu=px_rate,
                theta=px_r,
                zi_logits=px_dropout,
                scale=px_scale,
            )

        # we do not model the library size
        pl = None
        # prior
        pz = Normal(torch.zeros_like(z), torch.ones_like(z))

        return {
            MODULE_KEYS.PX_KEY: px,
            MODULE_KEYS.PL_KEY: pl,
            MODULE_KEYS.PZ_KEY: pz,
        }

    def loss(
        self,
        tensors: dict[str, torch.Tensor],
        inference_outputs: dict[str, torch.Tensor],
        generative_outputs: dict[str, torch.Tensor],
        lam_kl: torch.Tensor | float = 1.0,
        mode: int | None = None,
    ) -> LossOutput:
        x = tensors[REGISTRY_KEYS.X_KEY]
        n_obs = x.shape[0]
        n_var = x.shape[1]

        # data nll calculation
        reconst_loss = -generative_outputs[MODULE_KEYS.PX_KEY].log_prob(x).sum(-1)
        reconstruction_loss_norm = torch.mean(reconst_loss)

        # data kl div
        kl_divergence_z = kl_divergence(
            inference_outputs[MODULE_KEYS.QZ_KEY], generative_outputs[MODULE_KEYS.PZ_KEY]
        ).sum(dim=-1)

        kl_local_norm = torch.sum(kl_divergence_z) / (n_obs * n_var)

        loss = reconstruction_loss_norm + lam_kl * kl_local_norm

        ## graph inference
        mu_all = inference_outputs["mu_all"]
        logvar_all = inference_outputs["logvar_all"]
        v_all = inference_outputs["v_all"]

        return LossOutput(
            loss=loss,
            reconstruction_loss=reconst_loss,
            kl_local=kl_local_norm,
            extra_metrics={
                "mu_all": mu_all,
                "logvar_all": logvar_all,
                "v_all": v_all,
                "guidance_graph": self.guidance_graph,
            },
        )
