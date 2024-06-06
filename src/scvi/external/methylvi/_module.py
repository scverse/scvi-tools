"""PyTorch module for methylVI for single cell methylation data."""

from typing import Literal

import torch
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
        Human-readable list of methylation modalities (e.g. ["mCG", "mCH"])
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
        * ``'gene'`` - dispersion parameter of BetaBinomial is constant per gene across cells
        * ``'gene-cell'`` - dispersion can differ for every gene in every cell
    """

    def __init__(
        self,
        n_input: int,
        modalities: list[str],
        num_features_per_modality: list[int],
        n_batch: int = 0,
        n_hidden: int = 128,
        n_latent: int = 10,
        n_layers: int = 1,
        dropout_rate: float = 0.1,
        likelihood: Literal["betabinomial", "binomial"] = "betabinomial",
        dispersion: Literal["gene", "gene-cell"] = "gene",
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

        # z encoder goes from the n_input-dimensional data to an n_latent-d
        # latent space representation
        self.z_encoder = Encoder(
            n_input
            * 2,  # For each input region, we have both methylated counts and coverage --> x2
            n_latent,
            n_cat_list=cat_list,
            n_layers=n_layers,
            n_hidden=n_hidden,
            dropout_rate=dropout_rate,
            return_dist=True,
            var_activation=torch.nn.functional.softplus,  # Better numerical stability than exp
        )

        # decoder goes from n_latent-dimensional space to n_input-d data
        self.decoder = DecoderMETHYLVI(
            n_latent,
            n_input,
            n_cat_list=cat_list,
            n_layers=n_layers,
            n_hidden=n_hidden,
        )

        if self.dispersion == "gene":
            self.px_gamma = torch.nn.Parameter(torch.randn(n_input))
        elif self.dispersion == "gene-cell":
            pass

    def _get_methylation_features(self, tensors):
        if len(self.modalities) > 1:
            mc_keys = sorted(
                [x for x in tensors.keys() if x.endswith(f"_{METHYLVI_REGISTRY_KEYS.MC_KEY}")]
            )
            cov_keys = sorted(
                [x for x in tensors.keys() if x.endswith(f"_{METHYLVI_REGISTRY_KEYS.COV_KEY}")]
            )

            mc = torch.cat([tensors[x] for x in mc_keys], dim=1)
            cov = torch.cat([tensors[x] for x in cov_keys], dim=1)

        else:
            mc = tensors[METHYLVI_REGISTRY_KEYS.MC_KEY]
            cov = tensors[METHYLVI_REGISTRY_KEYS.COV_KEY]

        return mc, cov

    def _split_by_modality(self, tensor):
        tensor_by_modality = {}
        idxs = []
        for i in range(len(self.num_features_per_modality)):
            if i == 0:
                start = 0
            else:
                start = idxs[i - 1][1]
            end = start + self.num_features_per_modality[i]
            idxs.append((start, end))

        for i, modality in enumerate(self.modalities):
            tensor_by_modality[modality] = tensor[..., idxs[i][0] : idxs[i][1]]

        return tensor_by_modality

    def _get_methylation_features(self, tensors):
        if len(self.modalities) > 1:
            mc_keys = [f"{x}_{METHYLVI_REGISTRY_KEYS.MC_KEY}" for x in self.modalities]
            cov_keys = [f"{x}_{METHYLVI_REGISTRY_KEYS.COV_KEY}" for x in self.modalities]

            mc = torch.cat([tensors[x] for x in mc_keys], dim=1)
            cov = torch.cat([tensors[x] for x in cov_keys], dim=1)
        else:
            mc = tensors[METHYLVI_REGISTRY_KEYS.MC_KEY]
            cov = tensors[METHYLVI_REGISTRY_KEYS.COV_KEY]

        return mc, cov

    def _get_inference_input(self, tensors):
        """Parse the dictionary to get appropriate args"""
        mc, cov = self._get_methylation_features(tensors)

        batch_index = tensors[REGISTRY_KEYS.BATCH_KEY]

        input_dict = {
            METHYLVI_REGISTRY_KEYS.MC_KEY: mc,
            METHYLVI_REGISTRY_KEYS: cov,
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
        px_mu, px_gamma = self.decoder(z, batch_index)
        pz = Normal(torch.zeros_like(z), torch.ones_like(z))

        px_mu_by_modality = self._split_by_modality(px_mu)
        px_gamma_by_modality = self._split_by_modality(px_gamma)

        return {"px_mu": px_mu_by_modality, "px_gamma": px_gamma_by_modality, "pz": pz}

    def loss(
        self,
        tensors,
        inference_outputs,
        generative_outputs,
        kl_weight: float = 1.0,
    ):
        """Loss function."""
        mc, cov = self._get_methylation_features(tensors)

        px_mu = generative_outputs["px_mu"]
        px_gamma = generative_outputs["px_gamma"]

        px_mu = torch.concatenate([px_mu[modality] for modality in self.modalities], dim=1)
        px_gamma = torch.concatenate([px_gamma[modality] for modality in self.modalities], dim=1)

        if self.dispersion == "gene":
            px_gamma = torch.sigmoid(self.px_gamma)

        qz = inference_outputs["qz"]
        pz = generative_outputs["pz"]
        kl_divergence_z = kl(qz, pz).sum(dim=1)

        if self.likelihood == "binomial":
            dist = Binomial(probs=px_mu, total_count=cov)
        elif self.likelihood == "betabinomial":
            dist = BetaBinomial(mu=px_mu, gamma=px_gamma, total_count=cov)

        reconst_loss = -dist.log_prob(mc).sum(dim=-1)

        kl_local_for_warmup = kl_divergence_z

        weighted_kl_local = kl_weight * kl_local_for_warmup

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
    ) -> torch.Tensor:
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
            tensor with shape (n_cells, n_genes, n_samples)
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

        px_mu = generative_outputs["px_mu"]
        px_gamma = generative_outputs["px_gamma"]

        px_mu = torch.concatenate([px_mu[modality] for modality in self.modalities], dim=-1)
        px_gamma = torch.concatenate([px_gamma[modality] for modality in self.modalities], dim=-1)

        if self.dispersion == "gene":
            px_gamma = torch.sigmoid(self.px_gamma)

        mc, cov = self._get_methylation_features(tensors)

        if self.likelihood == "binomial":
            dist = Binomial(probs=px_mu, total_count=cov)
        elif self.likelihood == "betabinomial":
            dist = BetaBinomial(mu=px_mu, gamma=px_gamma, total_count=cov)

        if n_samples > 1:
            exprs = dist.sample()
            exprs = self._split_by_modality(exprs)
            for modality in self.modalities:
                exprs[modality] = (
                    exprs[modality].permute([1, 2, 0]).cpu()
                )  # Shape : (n_cells_batch, n_genes, n_samples)
        else:
            exprs = dist.sample().cpu()

        return exprs
