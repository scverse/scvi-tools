import torch
import torch.nn as nn
import torch.nn.functional as F
from torch.distributions import Normal
from torch.distributions import kl_divergence as kl

from scvi import REGISTRY_KEYS
from scvi.distributions import NegativeBinomial
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

        self.v = nn.Parameter(torch.randn(n_output, n_latent) * 0.01)

    def forward(
        self,
        u: torch.Tensor,
        l: torch.Tensor,
        batch_index: torch.Tensor,
    ) -> tuple[torch.Tensor, torch.Tensor]:
        # bring batch index in the right dimension
        if batch_index.dim() > 1:
            batch_index = batch_index.squeeze(-1)
        print("scale_shape: ", self.scale_lin.shape)
        print("bias_shape: ", self.bias.shape)
        # bring scale and bias in dimension [batch_size, n_genes] - row corresponds to batch index
        scale = F.softplus(self.scale_lin[batch_index])
        bias = self.bias[batch_index]
        log_theta = self.log_theta[batch_index]

        logit_mu = scale * (u @ self.v.T) + bias
        mu = F.softmax(logit_mu, dim=1) * l

        return mu, log_theta


class SPAGLUEVAE(BaseModuleClass):
    def __init__(
        self,
        n_inputs: list[int],
        n_batches: list[int],
        gene_likelihoods: list[str],
        n_latent_seq: int = 10,
        n_latent_spatial: int = 10,
        n_hidden: int = 128,
        n_layers: int = 2,
        dropout_rate: float = 0.1,
        **kwargs: dict,
    ) -> None:
        super().__init__()

        # Initialize attributes
        self.n_input_list = n_inputs
        self.n_batches_list = n_batches
        self.gene_likelihoods = gene_likelihoods

        # Define encoder networks
        self.z_encoder_diss = Encoder(
            n_input=n_inputs[0],
            n_output=n_latent_seq,
            n_hidden=n_hidden,
            n_layers=n_layers,
            dropout_rate=dropout_rate,
            **kwargs,
        )

        self.z_encoder_spa = Encoder(
            n_input=n_inputs[1],
            n_output=n_latent_spatial,
            n_hidden=n_hidden,
            n_layers=n_layers,
            dropout_rate=dropout_rate,
            **kwargs,
        )

        # Trainable linear decoders with batch variable
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

    def _get_inference_input(
        self, tensors: dict[str, torch.Tensor]
    ) -> dict[str, torch.Tensor | None]:
        """Get the input for the inference model."""
        return {
            "x": tensors[REGISTRY_KEYS.X_KEY],
        }

    def _get_generative_input(
        self, tensors: dict[str, torch.Tensor], inference_outputs: dict[str, torch.Tensor]
    ) -> dict[str, torch.Tensor]:
        """Get the input for the generative model."""
        z = inference_outputs["z"]
        library = inference_outputs["library"]
        batch_index = tensors[REGISTRY_KEYS.BATCH_KEY]
        y = tensors[REGISTRY_KEYS.LABELS_KEY]
        return {"z": z, "library": library, "batch_index": batch_index, "y": y}

    @auto_move_data
    def inference(
        self,
        x: torch.Tensor,
        mode: int | None = 0,
    ) -> dict[str, torch.Tensor]:
        """Run the inference model."""
        x_ = x

        # diss data
        if mode == 0:
            qm, qv, z = self.z_encoder_diss(x_)
        # spa data
        else:
            qm, qv, z = self.z_encoder_spa(x_)

        library = torch.sum(x, dim=1).view(-1, 1)

        # print("Inference Output")
        # print("z: ", z.shape)
        # print("qm: ", qm.shape)
        # print("qv: ", qv.shape)
        # print("library: ", library.shape)

        return {"qm": qm, "qv": qv, "z": z, "library": library}

    @auto_move_data
    def generative(
        self,
        z: torch.Tensor,
        library: torch.Tensor,
        batch_index: torch.Tensor | None = None,
        y: torch.Tensor | None = None,
        mode: int | None = 0,
    ) -> dict[str, torch.Tensor]:
        """Run the generative model."""
        # diss data
        if mode == 0:
            mu, log_theta = self.z_decoder_diss(z, library, batch_index)
        # spa data
        else:
            mu, log_theta = self.z_decoder_spa(z, library, batch_index)

        # print("Generative Output: ")
        # print("mu: ", mu.shape)
        # print("log_theta: ", log_theta.shape)

        return {
            "mu": mu,
            "log_theta": log_theta,
        }

    def reconstruction_loss(
        self,
        x: torch.Tensor,
        mu: torch.Tensor,
        log_theta: torch.Tensor,
        mode: int,
    ) -> torch.Tensor:
        reconstruction_loss = None

        if self.gene_likelihoods[mode] == "nb":
            reconstruction_loss = (
                -NegativeBinomial(log_theta.exp(), logits=(mu + EPS).log() - log_theta)
                .log_prob(x)
                .sum(dim=-1)
            )

        return reconstruction_loss

    def loss(
        self,
        tensors: dict[str, torch.Tensor],
        inference_outputs: dict[str, torch.Tensor],
        generative_outputs: dict[str, torch.Tensor],
        mode: int | None = None,
    ) -> LossOutput:
        x = tensors[REGISTRY_KEYS.X_KEY]

        # reconstruction loss
        mu = generative_outputs["mu"]
        log_theta = generative_outputs["log_theta"]

        reconstruction_loss = self.reconstruction_loss(
            x,
            mu,
            log_theta,
            mode,
        )

        # KL Divergence
        qm = inference_outputs["qm"]
        qv = inference_outputs["qv"]

        qz = Normal(loc=qm, scale=torch.sqrt(qv))
        pz = Normal(loc=torch.zeros_like(qm), scale=torch.ones_like(qv))

        kl_local = kl(qz, pz).sum(dim=1)

        n_obs = x.shape[0]
        n_var = x.shape[1]

        # normalization not with latent dimension?
        # n_var_latent = qm.shape[1]

        reconstruction_loss_norm = torch.mean(reconstruction_loss)
        kl_local_norm = torch.sum(kl_local) / (n_obs * n_var)

        # print("reconstruction loss: ", reconstruction_loss_norm)
        # print("kl local: ", kl_local_norm)

        loss = reconstruction_loss_norm + kl_local_norm

        print("Total Loss: ", loss)

        return LossOutput(loss=loss, reconstruction_loss=reconstruction_loss, kl_local=kl_local)
