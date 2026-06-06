from __future__ import annotations

import torch
import torch.nn as nn
import torch.nn.functional as F
from torch.distributions import Normal, kl_divergence

from scvi.module.base import BaseModuleClass, LossOutput, auto_move_data
from scvi.nn import DecoderSCVI, Encoder


class SCPOLIModule(BaseModuleClass):
    """
    Module for scPoli.

    Parameters
    ----------
    n_input
        Number of input genes.
    n_conditions
        Number of unique samples/batches.
    n_cell_types
        Number of cell types (for prototype layer).
    n_latent
        Dimensionality of latent space.
    embedding_dim
        Dimensionality of sample embeddings.
    hidden_dims
        Hidden layer sizes for encoder/decoder.
    recon_loss
        Reconstruction loss: "nb" or "zinb".
    """

    def __init__(
        self,
        n_input: int,
        n_conditions: int,
        n_cell_types: int,
        n_latent: int = 10,
        embedding_dim: int = 5,
        hidden_dims: list[int] | None = None,
        recon_loss: str = "nb",
    ):
        super().__init__()

        if hidden_dims is None:
            hidden_dims = [256, 64]

        self.n_latent = n_latent
        self.n_cell_types = n_cell_types
        self.recon_loss = recon_loss
        self.embed = nn.Embedding(n_conditions, embedding_dim)
        self.encoder = Encoder(
            n_input=n_input + embedding_dim,
            n_output=n_latent,
            n_hidden=hidden_dims[0],
            n_layers=len(hidden_dims),
            dropout_rate=0.1,
        )
        self.decoder = DecoderSCVI(
            n_input=n_latent + embedding_dim,
            n_output=n_input,
            n_hidden=hidden_dims[0],
            n_layers=len(hidden_dims),
        )
        self.prototypes = nn.Parameter(torch.randn(n_cell_types, n_latent))

        @auto_move_data
        def inference(self, x, condition_id, **kwargs):
            """Encoder forward pass."""
            # Get sample embedding for each cell
            cond_emb = self.embed(condition_id.squeeze(-1))

            # Concatenate gene expression with sample embedding
            x_cond = torch.cat([x, cond_emb], dim=-1)

            # Encode to latent space
            qz_m, qz_v, z = self.encoder(x_cond)

            return {z: z, qz_m: qz_m, qz_v: qz_v, cond_emb: cond_emb}

        @auto_move_data
        def generative(self, z, cond_emb, **kwargs):
            """Decoder forward pass."""
            # Concatenate latent with sample embedding
            z_cond = torch.cat([z, cond_emb], dim=-1)

            # Decode back to gene expression
            px_scale, px_r, px_rate, px_dropout = self.decoder(
                self.recon_loss,
                z_cond,
                torch.ones(z.shape[0], 1).to(z.device),  # library size
            )
            return {
                px_scale: px_scale,
                px_r: px_r,
                px_rate: px_rate,
                px_dropout: px_dropout,
            }

        def loss(
            self,
            tensors,
            inference_outputs,
            generative_outputs,
            kl_weight: float = 1.0,
            proto_weight: float = 1.0,
        ):
            """Compute ELBO + prototype classification loss."""
            from scvi.distributions._negative_binomial import (
                NegativeBinomial,
                ZeroInflatedNegativeBinomial,
            )

            x = tensors["X"]
            labels = tensors.get("labels", None)

            z = inference_outputs["z"]
            qz_m = inference_outputs["qz_m"]
            qz_v = inference_outputs["qz_v"]

            # 1. Reconstruction loss
            if self.recon_loss == "zinb":
                recon = (
                    ZeroInflatedNegativeBinomial(
                        mu=generative_outputs["px_rate"],
                        theta=generative_outputs["px_r"],
                        zi_logits=generative_outputs["px_dropout"],
                    )
                    .log_prob(x)
                    .sum(-1)
                )
            else:
                recon = (
                    NegativeBinomial(
                        mu=generative_outputs["px_rate"],
                        theta=generative_outputs["px_r"],
                    )
                    .log_prob(x)
                    .sum(-1)
                )

            # 2. KL divergence
            kl = kl_divergence(
                Normal(qz_m, qz_v.sqrt()),
                Normal(torch.zeros_like(qz_m), torch.ones_like(qz_v)),
            ).sum(-1)
            # 3. Prototype loss (only when labels are available)
            proto_loss = torch.tensor(0.0, device=z.device)
            if labels is not None:
                # Distance from each cell to every prototype
                dists = torch.cdist(z, self.prototypes)  # (n_cells, n_cell_types)
                # Cell should be closest to its own cell type prototype
                proto_loss = F.cross_entropy(-dists, labels.squeeze(-1).long())

                # Total loss
            elbo = -recon + kl_weight * kl
            loss = elbo.mean() + proto_weight * proto_loss

            return LossOutput(
                loss=loss,
                reconstruction_loss=-recon,
                kl_local=kl,
            )
