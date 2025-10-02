# -*- coding: utf-8 -*-
"""Main module."""
from typing import Callable, Iterable, Optional, Literal

import numpy as np
import torch
import torch.nn.functional as F
from torch.nn.functional import one_hot
from torch.distributions import Normal, Poisson

from scvi import REGISTRY_KEYS
from scvi.distributions import NegativeBinomial, ZeroInflatedNegativeBinomial
from scvi.module.base import BaseModuleClass, LossOutput, auto_move_data
from scvi.nn import DecoderSCVI, Encoder
from scvi.module import VAE

torch.backends.cudnn.benchmark = True


# WAE model
class WAE(VAE):
    """
    Wasserstein auto-encoder model.

    This is an implementation of the WAE model described in [Ergen22]_

    Parameters
    ----------
    n_input
        Number of input genes
    n_batch
        Number of batches, if 0, no batch correction is performed.
    n_labels
        Number of labels
    n_hidden
        Number of nodes per hidden layer
    n_latent
        Dimensionality of the latent space
    n_layers
        Number of hidden layers used for encoder and decoder NNs
    n_continuous_cov
        Number of continuous covariates
    n_cats_per_cov
        Number of categories for each extra categorical covariate
    dropout_rate
        Dropout rate for neural networks
    dispersion
        One of the following

        * ``'gene'`` - dispersion parameter of NB is constant per gene across cells
        * ``'gene-batch'`` - dispersion can differ between different batches
        * ``'gene-label'`` - dispersion can differ between different labels
        * ``'gene-cell'`` - dispersion can differ for every gene in every cell
    log_variational
        Log(data+1) prior to encoding for numerical stability. Not normalization.
    gene_likelihood
        One of

        * ``'nb'`` - Negative binomial distribution
        * ``'zinb'`` - Zero-inflated negative binomial distribution
        * ``'poisson'`` - Poisson distribution
    latent_distribution
        One of

        * ``'normal'`` - Isotropic normal
        * ``'ln'`` - Logistic normal with normal params N(0, 1)
    encode_covariates
        Whether to concatenate covariates to expression in encoder
    deeply_inject_covariates
        Whether to concatenate covariates into output of hidden layers in encoder/decoder. This option
        only applies when `n_layers` > 1. The covariates are concatenated to the input of subsequent hidden layers.
    use_layer_norm
        Whether to use layer norm in layers
    use_size_factor_key
        Use size_factor AnnDataField defined by the user as scaling factor in mean of conditional distribution.
        Takes priority over `use_observed_lib_size`.
    var_activation
        Callable used to ensure positivity of the variational distributions' variance.
        When `None`, defaults to `torch.exp`.

    """

    def __init__(
        self,
        n_input: int,
        n_batch: int = 0,
        n_labels: int = 1,
        n_hidden: int = 128,
        n_latent: int = 48,
        n_layers: int = 2,
        n_continuous_cov: int = 0,
        n_cats_per_cov: Optional[Iterable[int]] = None,
        dropout_rate: float = 0.1,
        dispersion: str = "gene",
        log_variational: bool = True,
        gene_likelihood: str = "nb",
        latent_distribution: str = "normal",
        encode_covariates: bool = False,
        deeply_inject_covariates: bool = True,
        use_batch_norm: Literal["encoder", "decoder", "none", "both"] = "none",
        use_layer_norm: Literal["encoder", "decoder", "none", "both"] = "both",
        use_size_factor_key: bool = False,
        var_activation: Optional[Callable] = torch.sigmoid,
        num_sample: int = 8,
        weighting_genes: Optional[np.ndarray] = None,
    ):
        super().__init__(
            n_input=n_input,
            n_batch=n_batch,
            n_labels=n_labels,
            n_hidden=n_hidden,
            n_latent=n_latent,
            n_layers=n_layers,
            n_continuous_cov=n_continuous_cov,
            n_cats_per_cov=n_cats_per_cov,
            dropout_rate=dropout_rate,
            dispersion=dispersion,
            log_variational=log_variational,
            gene_likelihood=gene_likelihood,
            latent_distribution=latent_distribution,
            encode_covariates=encode_covariates,
            deeply_inject_covariates=deeply_inject_covariates,
            use_batch_norm=use_batch_norm,
            use_layer_norm=use_layer_norm,
            use_size_factor_key=use_size_factor_key,
            var_activation=var_activation,
        )
        self.num_sample = num_sample

        self.register_buffer(
            "weighting_genes",
            torch.Tensor(weighting_genes) if weighting_genes is not None else torch.ones(n_input)
        )

    def loss(
        self,
        tensors,
        inference_outputs,
        generative_outputs,
        n_obs,
        weighting_random_encoder=0,
        wae_lambda=0,
        kl_weight: float = 1.0,  # Downscales mmd loss for the first epochs
    ):
        x = tensors[REGISTRY_KEYS.X_KEY]
        qz = inference_outputs[MODULE_KEYS.QZ_KEY]
        pz = generative_outputs[MODULE_KEYS.PZ_KEY]
        reconst_loss = - (self.weighting_genes * generative_outputs["px"].log_prob(x)).sum(-1)

        randomization_loss = torch.abs(torch.log(qz.var)).sum(dim=-1)
        mmd_loss = self.get_mmd_loss(qz, pz, n_obs)
        kl_divergence_l = 0.0

        loss = torch.mean(
            reconst_loss
            + kl_weight * (
                weighting_random_encoder * torch.mean(randomization_loss)
                + wae_lambda * mmd_loss
            )
            + kl_divergence_l
        )

        return LossOutput(
            loss=loss,
            reconstruction_loss=reconst_loss,
            kl_local=wae_lambda*mmd_loss + weighting_random_encoder*randomization_loss,
        )

    def get_mmd_loss(self, qz, pz, total_number_cells) -> torch.Tensor:
        batch_size = pz.shape[0] // self.num_sample
        latent_size = qz.shape[1]
        z_mmd = qz.rsample([self.num_sample])
        p_z = pz.rsample()

        distances_pz = torch.square(torch.cdist(p_z, p_z, p=2))
        distances_qz = torch.square(torch.cdist(z_mmd, z_mmd, p=2))
        distances_pz_qz = torch.square(torch.cdist(p_z, z_mmd, p=2))

        def _block_diagonal() -> torch.Tensor:
            first_idx = torch.repeat_interleave(
                torch.arange(0, self.num_sample * batch_size), repeats=self.num_sample
            ).to(qz.device)

            second_idx = (
                torch.reshape(
                    torch.arange(0, self.num_sample * batch_size), (-1, self.num_sample)
                )
                .repeat(1, self.num_sample)
                .flatten()
                .to(qz.device)
            )

            idx = torch.stack((first_idx, second_idx), dim=1)

            block = (
                torch.sparse_coo_tensor(
                    idx.t(),
                    torch.ones_like(idx[:, 0]),
                    torch.Size(
                        [batch_size * self.num_sample, batch_size * self.num_sample]
                    ),
                )
                .to_dense()
                .to(qz.device)
            )

            return block

        mask = _block_diagonal().to(qz.device)

        eye_matrix = torch.subtract(torch.ones(1), torch.eye(batch_size)).to(qz.device)
        second_eye = torch.subtract(
            mask, torch.eye((batch_size * self.num_sample)).to(qz.device)
        )

        cbase = 2.0 * latent_size
        stat = 0.0

        for scale in [0.1, 0.5, 1.0, 10.0]:
            cval = cbase * scale
            res1 = cval / (cval + distances_pz)
            res1 = torch.multiply(res1, eye_matrix)
            res1 = torch.sum(res1) / (batch_size - 1)

            res2 = cval / (cval + distances_pz_qz)
            res2 = torch.sum(res2) / (batch_size * self.num_sample)

            res3 = cval / (cval + distances_qz)
            res3 = torch.sum(torch.multiply(res3, 1.0 - mask))
            res3 = res3 / (batch_size - 1) / (self.num_sample**2)

            res4 = cval / (cval + distances_qz)
            res4 = torch.multiply(res4, second_eye)
            res4 = (
                torch.sum(res4)
                / self.num_sample
                / (self.num_sample - 1.0)
                * (batch_size / total_number_cells)
            )

            stat += res1 - 2 * res2 + res3 + res4

        return stat
