import torch
from torch.distributions import Normal, kl_divergence

from scvi import REGISTRY_KEYS
from scvi.distributions import NegativeBinomial, ZeroInflatedNegativeBinomial
from scvi.external.spaglue import GraphEncoder_glue, NBDataDecoderWB
from scvi.module._constants import MODULE_KEYS
from scvi.module.base import BaseModuleClass, LossOutput, auto_move_data
from scvi.nn import Encoder


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
            n_input=n_inputs["diss"],
            n_output=n_latent_seq,
            n_hidden=n_hidden,
            n_layers=n_layers,
            dropout_rate=dropout_rate,
            return_dist=True,
            # **kwargs,
        )

        self.z_encoder_spa = Encoder(
            n_input=n_inputs["spatial"],
            n_output=n_latent_spatial,
            n_hidden=n_hidden,
            n_layers=n_layers,
            dropout_rate=dropout_rate,
            return_dist=True,
            # **kwargs,
        )

        self.z_decoder_diss = NBDataDecoderWB(
            n_output=n_inputs["diss"],
            n_latent=n_latent_seq,
            n_batches=n_batches["diss"],
        )

        self.z_decoder_spa = NBDataDecoderWB(
            n_output=n_inputs["spatial"],
            n_latent=n_latent_spatial,
            n_batches=n_batches["spatial"],
        )

        self.graph_encoder = GraphEncoder_glue(
            vnum=n_inputs["diss"] + n_inputs["spatial"],
            out_features=50,
        )

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
        mode: str | None = None,
    ) -> dict[str, torch.Tensor]:
        x_ = x
        library = torch.log(x.sum(1)).unsqueeze(1)

        graph = self.guidance_graph
        device = x.device
        graph = graph.to(device)

        # whole embedding is calculated
        v_all, mu_all, logvar_all = self.graph_encoder(graph.edge_index)

        # embedding for modality is extracted to be used for decoder input
        if mode == "diss":
            v = v_all[getattr(graph, f"{mode}_indices")]
            other_mode = [m for m in ["diss", "spatial"] if m != mode][0]
            v_other_mod = v_all[getattr(graph, f"{other_mode}_indices")]
        elif mode == "spatial":
            v = v_all[getattr(graph, f"{mode}_indices")]
            other_mode = [m for m in ["diss", "spatial"] if m != mode][0]
            v_other_mod = v_all[getattr(graph, f"{other_mode}_indices")]
        else:
            raise ValueError("Invalid mode: must be diss or spatial.")

        # diss data
        if mode == "diss":
            qz, z = self.z_encoder_diss(x_)
        # spa data
        else:
            qz, z = self.z_encoder_spa(x_)

        return {
            MODULE_KEYS.QZ_KEY: qz,
            MODULE_KEYS.Z_KEY: z,
            MODULE_KEYS.LIBRARY_KEY: library,
            "v": v,
            "v_other": v_other_mod,
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
        mode: str | None = 0,
    ) -> dict[str, torch.Tensor]:
        """Run the generative model."""
        EPS = 1e-8

        # diss data
        if mode == "diss":
            px_scale, px_r, px_rate, px_dropout = self.z_decoder_diss(z, library, batch_index, v)

        # spa data
        elif mode == "spatial":
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
        # pl = None
        # prior
        pz = Normal(torch.zeros_like(z), torch.ones_like(z))

        return {
            MODULE_KEYS.PX_KEY: px,
            # MODULE_KEYS.PL_KEY: pl,
            MODULE_KEYS.PZ_KEY: pz,
            "px_rate": px_rate,
        }

    def loss(
        self,
        tensors: dict[str, torch.Tensor],
        inference_outputs: dict[str, torch.Tensor],
        generative_outputs: dict[str, torch.Tensor],
        lam_kl: torch.Tensor | float = 1.0,
        lam_data: torch.Tensor | float = 1.0,
        mode: str | None = None,
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

        loss = lam_data * reconstruction_loss_norm + lam_kl * kl_local_norm

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
