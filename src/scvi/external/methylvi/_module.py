"""PyTorch module for methylVI for single cell methylation data."""

from collections.abc import Iterable
from typing import Literal

import torch
import torch.nn as nn
from torch.distributions import Binomial, Normal
from torch.distributions import kl_divergence as kl

from scvi import REGISTRY_KEYS
from scvi.distributions import BetaBinomial
from scvi.external.methylvi import METHYLVI_REGISTRY_KEYS, DecoderMETHYLVI
from scvi.module.base import BaseModuleClass, LossOutput, auto_move_data
from scvi.nn import Encoder

TensorDict = dict[str, torch.Tensor]


class METHYLVAE(BaseModuleClass):
    """
    Pytorch module for methylVI.

    Parameters
    ----------
    n_input
        Total number of input genomic regions
    modalities
        List of methylation modalities (e.g. ["mCG", "mCH"])
    num_features_per_modality
        Number of features corresponding to each modality
    n_batch
        Number of batches, if 0, no batch correction is performed.
    n_hidden
        Number of nodes per hidden layer
    n_latent
        Dimensionality of the latent space
    n_layers
        Number of hidden layers used for encoder and decoder NNs
    dropout_rate
        Dropout rate for neural networks
    likelihood
        One of
        * ``'betabinomial'`` - BetaBinomial distribution
        * ``'binomial'`` - Binomial distribution
    dispersion
        One of the following
        * ``'region'`` - dispersion parameter of BetaBinomial is constant per region across cells
        * ``'region-cell'`` - dispersion can differ for every regionin every cell
    """

    def __init__(
        self,
        n_input: int,
        modalities: Iterable[str],
        num_features_per_modality: Iterable[int],
        n_batch: int = 0,
        n_hidden: int = 128,
        n_latent: int = 10,
        n_layers: int = 1,
        dropout_rate: float = 0.1,
        likelihood: Literal["betabinomial", "binomial"] = "betabinomial",
        dispersion: Literal["region", "region-cell"] = "region",
    ):
        super().__init__()
        self.n_latent = n_latent
        self.n_batch = n_batch

        self.latent_distribution = "normal"
        self.dispersion = dispersion
        self.likelihood = likelihood
        self.modalities = modalities
        self.num_features_per_modality = num_features_per_modality

        cat_list = [n_batch]

        self.z_encoder = Encoder(
            n_input * 2,  # Methylated counts and coverage for each feature --> x2
            n_latent,
            n_cat_list=cat_list,
            n_layers=n_layers,
            n_hidden=n_hidden,
            dropout_rate=dropout_rate,
            return_dist=True,
            var_activation=torch.nn.functional.softplus,  # Better numerical stability than exp
        )

        self.decoders = nn.ModuleDict()
        for modality, num_features in zip(modalities, num_features_per_modality):
            self.decoders[modality] = DecoderMETHYLVI(
                n_latent,
                num_features,
                n_cat_list=cat_list,
                n_layers=n_layers,
                n_hidden=n_hidden,
            )

        if self.dispersion == "region":
            self.px_gamma = torch.nn.ParameterDict(
                {
                    modality: nn.Parameter(torch.randn(num_features))
                    for (modality, num_features) in zip(modalities, num_features_per_modality)
                }
            )

    def _get_inference_input(self, tensors):
        """Parse the dictionary to get appropriate args"""
        mc = torch.cat(
            [tensors[f"{x}_{METHYLVI_REGISTRY_KEYS.MC_KEY}"] for x in sorted(self.modalities)],
            dim=1,
        )
        cov = torch.cat(
            [tensors[f"{x}_{METHYLVI_REGISTRY_KEYS.COV_KEY}"] for x in sorted(self.modalities)],
            dim=1,
        )

        batch_index = tensors[REGISTRY_KEYS.BATCH_KEY]

        input_dict = {
            METHYLVI_REGISTRY_KEYS.MC_KEY: mc,
            METHYLVI_REGISTRY_KEYS.COV_KEY: cov,
            "batch_index": batch_index,
        }
        return input_dict

    def _get_generative_input(self, tensors, inference_outputs):
        z = inference_outputs["z"]
        batch_index = tensors[REGISTRY_KEYS.BATCH_KEY]

        input_dict = {
            "z": z,
            "batch_index": batch_index,
        }
        return input_dict

    @auto_move_data
    def inference(self, mc, cov, batch_index, n_samples=1):
        """
        High level inference method.

        Runs the inference (encoder) model.
        """
        # log the inputs to the variational distribution for numerical stability
        mc_ = torch.log(1 + mc)
        cov_ = torch.log(1 + cov)

        # get variational parameters via the encoder networks
        # we input both the methylated reads (mc) and coverage (cov)
        total_input = torch.cat((mc_, cov_), dim=-1)
        qz, z = self.z_encoder(total_input, batch_index)
        if n_samples > 1:
            z = qz.sample((n_samples,))

        outputs = {"z": z, "qz": qz}
        return outputs

    @auto_move_data
    def generative(self, z, batch_index):
        """Runs the generative model."""
        # form the parameters of the BetaBinomial likelihood
        px_mu, px_gamma = {}, {}
        for modality in self.modalities:
            px_mu[modality], px_gamma[modality] = self.decoders[modality](
                self.dispersion, z, batch_index
            )

        pz = Normal(torch.zeros_like(z), torch.ones_like(z))

        return {"px_mu": px_mu, "px_gamma": px_gamma, "pz": pz}

    def loss(
        self,
        tensors,
        inference_outputs,
        generative_outputs,
        kl_weight: float = 1.0,
    ):
        """Loss function."""
        qz = inference_outputs["qz"]
        pz = generative_outputs["pz"]
        kl_divergence_z = kl(qz, pz).sum(dim=1)

        kl_local_for_warmup = kl_divergence_z

        weighted_kl_local = kl_weight * kl_local_for_warmup

        minibatch_size = qz.loc.size()[0]
        reconst_loss = torch.zeros(minibatch_size).to(self.device)

        for modality in self.modalities:
            px_mu = generative_outputs["px_mu"][modality]
            px_gamma = generative_outputs["px_gamma"][modality]
            mc = tensors[f"{modality}_{METHYLVI_REGISTRY_KEYS.MC_KEY}"]
            cov = tensors[f"{modality}_{METHYLVI_REGISTRY_KEYS.COV_KEY}"]

            if self.dispersion == "region":
                px_gamma = torch.sigmoid(self.px_gamma[modality])

            if self.likelihood == "binomial":
                dist = Binomial(probs=px_mu, total_count=cov)
            elif self.likelihood == "betabinomial":
                dist = BetaBinomial(mu=px_mu, gamma=px_gamma, total_count=cov)

            reconst_loss += -dist.log_prob(mc).sum(dim=-1)

        loss = torch.mean(reconst_loss + weighted_kl_local)

        kl_local = {"kl_divergence_z": kl_divergence_z}
        return LossOutput(
            loss=loss,
            reconstruction_loss=reconst_loss,
            kl_local=kl_local,
        )

    @torch.no_grad()
    def sample(
        self,
        tensors,
        n_samples=1,
    ) -> dict[torch.Tensor]:
        r"""
        Generate observation samples from the posterior predictive distribution.

        The posterior predictive distribution is written as :math:`p(\hat{x} \mid x)`.

        Parameters
        ----------
        tensors
            Tensors dict
        n_samples
            Number of required samples for each cell

        Returns
        -------
        x_new
            tensor with shape (n_cells, n_regions, n_samples)
        """
        inference_kwargs = {"n_samples": n_samples}
        (
            _,
            generative_outputs,
        ) = self.forward(
            tensors,
            inference_kwargs=inference_kwargs,
            compute_loss=False,
        )

        exprs = {}
        for modality in self.modalities:
            px_mu = generative_outputs["px_mu"][modality]
            px_gamma = generative_outputs["px_gamma"][modality]
            cov = tensors[f"{modality}_{METHYLVI_REGISTRY_KEYS.COV_KEY}"]

            if self.dispersion == "region":
                px_gamma = torch.sigmoid(self.px_gamma[modality])

            if self.likelihood == "binomial":
                dist = Binomial(probs=px_mu, total_count=cov)
            elif self.likelihood == "betabinomial":
                dist = BetaBinomial(mu=px_mu, gamma=px_gamma, total_count=cov)

            if n_samples > 1:
                exprs_ = dist.sample()
                exprs[modality] = exprs_.permute(
                    [1, 2, 0]
                ).cpu()  # Shape : (n_cells_batch, n_regions, n_samples)
            else:
                exprs[modality] = dist.sample().cpu()

        return exprs
