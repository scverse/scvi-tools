# -*- coding: utf-8 -*-
"""Main module."""
from typing import Dict, Optional, Tuple, Union

import numpy as np
import torch
import torch.nn as nn
import torch.nn.functional as F
from torch.distributions import Normal
from torch.distributions import kl_divergence as kl

from scvi._compat import Literal
from scvi.core.distributions import (
    NegativeBinomial,
    NegativeBinomialMixture,
    ZeroInflatedNegativeBinomial,
)

from ._base import DecoderTOTALVI, EncoderTOTALVI
from .utils import one_hot

torch.backends.cudnn.benchmark = True


# VAE model
class TOTALVAE(nn.Module):
    """
    Total variational inference for CITE-seq data.

    Implements the totalVI model of [GayosoSteier20]_.

    Parameters
    ----------
    n_input_genes
        Number of input genes
    n_input_proteins
        Number of input proteins
    n_batch
        Number of batches
    n_labels
        Number of labels
    n_hidden
        Number of nodes per hidden layer for encoder and decoder
    n_latent
        Dimensionality of the latent space
    n_layers
        Number of hidden layers used for encoder and decoder NNs
    dropout_rate
        Dropout rate for neural networks
    gene_dispersion
        One of the following

        * ``'gene'`` - genes_dispersion parameter of NB is constant per gene across cells
        * ``'gene-batch'`` - genes_dispersion can differ between different batches
        * ``'gene-label'`` - genes_dispersion can differ between different labels
    protein_dispersion
        One of the following

        * ``'protein'`` - protein_dispersion parameter is constant per protein across cells
        * ``'protein-batch'`` - protein_dispersion can differ between different batches NOT TESTED
        * ``'protein-label'`` - protein_dispersion can differ between different labels NOT TESTED
    log_variational
        Log(data+1) prior to encoding for numerical stability. Not normalization.
    gene_likelihood
        One of

        * ``'nb'`` - Negative binomial distribution
        * ``'zinb'`` - Zero-inflated negative binomial distribution
    latent_distribution
        One of

        * ``'normal'`` - Isotropic normal
        * ``'ln'`` - Logistic normal with normal params N(0, 1)
    protein_batch_mask
        Dictionary where each key is a batch code, and value is for each protein, whether it was observed or not.
    encode_covariates
        Whether to concatenate covariates to expression in encoder
    protein_background_prior_mean
        Array of proteins by batches, the prior initialization for the protein background mean (log scale)
    protein_background_prior_scale
        Array of proteins by batches, the prior initialization for the protein background scale (log scale)
    use_observed_lib_size
        Use observed library size for RNA as scaling factor in mean of conditional distribution
    """

    def __init__(
        self,
        n_input_genes: int,
        n_input_proteins: int,
        n_batch: int = 0,
        n_labels: int = 0,
        n_hidden: int = 256,
        n_latent: int = 20,
        n_layers_encoder: int = 2,
        n_layers_decoder: int = 1,
        dropout_rate_decoder: float = 0.2,
        dropout_rate_encoder: float = 0.2,
        gene_dispersion: str = "gene",
        protein_dispersion: str = "protein",
        log_variational: bool = True,
        gene_likelihood: str = "nb",
        latent_distribution: str = "normal",
        protein_batch_mask: Dict[Union[str, int], np.ndarray] = None,
        encode_covariates: bool = True,
        protein_background_prior_mean: Optional[np.ndarray] = None,
        protein_background_prior_scale: Optional[np.ndarray] = None,
        use_observed_lib_size: bool = True,
        use_batch_norm: Literal["encoder", "decoder", "none", "both"] = "both",
        use_layer_norm: Literal["encoder", "decoder", "none", "both"] = "none",
    ):
        super().__init__()
        self.gene_dispersion = gene_dispersion
        self.n_latent = n_latent
        self.log_variational = log_variational
        self.gene_likelihood = gene_likelihood
        self.n_batch = n_batch
        self.n_labels = n_labels
        self.n_input_genes = n_input_genes
        self.n_input_proteins = n_input_proteins
        self.protein_dispersion = protein_dispersion
        self.latent_distribution = latent_distribution
        self.protein_batch_mask = protein_batch_mask
        self.use_observed_lib_size = use_observed_lib_size

        # parameters for prior on rate_back (background protein mean)
        if protein_background_prior_mean is None:
            if n_batch > 0:
                self.background_pro_alpha = torch.nn.Parameter(
                    torch.randn(n_input_proteins, n_batch)
                )
                self.background_pro_log_beta = torch.nn.Parameter(
                    torch.clamp(torch.randn(n_input_proteins, n_batch), -10, 1)
                )
            else:
                self.background_pro_alpha = torch.nn.Parameter(
                    torch.randn(n_input_proteins)
                )
                self.background_pro_log_beta = torch.nn.Parameter(
                    torch.clamp(torch.randn(n_input_proteins), -10, 1)
                )
        else:
            if protein_background_prior_mean.shape[1] == 1 and n_batch != 1:
                init_mean = protein_background_prior_mean.ravel()
                init_scale = protein_background_prior_scale.ravel()
            else:
                init_mean = protein_background_prior_mean
                init_scale = protein_background_prior_scale
            self.background_pro_alpha = torch.nn.Parameter(
                torch.from_numpy(init_mean.astype(np.float32))
            )
            self.background_pro_log_beta = torch.nn.Parameter(
                torch.log(torch.from_numpy(init_scale.astype(np.float32)))
            )

        if self.gene_dispersion == "gene":
            self.px_r = torch.nn.Parameter(torch.randn(n_input_genes))
        elif self.gene_dispersion == "gene-batch":
            self.px_r = torch.nn.Parameter(torch.randn(n_input_genes, n_batch))
        elif self.gene_dispersion == "gene-label":
            self.px_r = torch.nn.Parameter(torch.randn(n_input_genes, n_labels))
        else:  # gene-cell
            pass

        if self.protein_dispersion == "protein":
            self.py_r = torch.nn.Parameter(2 * torch.rand(self.n_input_proteins))
        elif self.protein_dispersion == "protein-batch":
            self.py_r = torch.nn.Parameter(
                2 * torch.rand(self.n_input_proteins, n_batch)
            )
        elif self.protein_dispersion == "protein-label":
            self.py_r = torch.nn.Parameter(
                2 * torch.rand(self.n_input_proteins, n_labels)
            )
        else:  # protein-cell
            pass

        use_batch_norm_encoder = use_batch_norm == "encoder" or use_batch_norm == "both"
        use_batch_norm_decoder = use_batch_norm == "decoder" or use_batch_norm == "both"
        use_layer_norm_encoder = use_layer_norm == "encoder" or use_layer_norm == "both"
        use_layer_norm_decoder = use_layer_norm == "decoder" or use_layer_norm == "both"

        # z encoder goes from the n_input-dimensional data to an n_latent-d
        # latent space representation
        self.encoder = EncoderTOTALVI(
            n_input_genes + self.n_input_proteins,
            n_latent,
            n_layers=n_layers_encoder,
            n_cat_list=[n_batch] if encode_covariates else None,
            n_hidden=n_hidden,
            dropout_rate=dropout_rate_encoder,
            distribution=latent_distribution,
            use_batch_norm=use_batch_norm_encoder,
            use_layer_norm=use_layer_norm_encoder,
        )
        self.decoder = DecoderTOTALVI(
            n_latent,
            n_input_genes,
            self.n_input_proteins,
            n_layers=n_layers_decoder,
            n_cat_list=[n_batch],
            n_hidden=n_hidden,
            dropout_rate=dropout_rate_decoder,
            use_batch_norm=use_batch_norm_decoder,
            use_layer_norm=use_layer_norm_decoder,
        )

    def sample_from_posterior_z(
        self,
        x: torch.Tensor,
        y: torch.Tensor,
        batch_index: Optional[torch.Tensor] = None,
        give_mean: bool = False,
        n_samples: int = 5000,
    ) -> torch.Tensor:
        """
        Access the tensor of latent values from the posterior.

        Parameters
        ----------
        x
            tensor of values with shape ``(batch_size, n_input_genes)``
        y
            tensor of values with shape ``(batch_size, n_input_proteins)``
        batch_index
            tensor of batch indices
        give_mean
            Whether to sample, or give mean of distribution
        n_samples
            Number of samples for monte carlo estimation

        Returns
        -------
        type
            tensor of shape ``(batch_size, n_latent)``
        """
        if self.log_variational:
            x = torch.log(1 + x)
            y = torch.log(1 + y)
        qz_m, qz_v, _, _, latent, _ = self.encoder(
            torch.cat((x, y), dim=-1), batch_index
        )
        z = latent["z"]
        if give_mean:
            if self.latent_distribution == "ln":
                samples = Normal(qz_m, qz_v.sqrt()).sample([n_samples])
                z = self.encoder.z_transformation(samples)
                z = z.mean(dim=0)
            else:
                z = qz_m
        return z

    def sample_from_posterior_l(
        self,
        x: torch.Tensor,
        y: torch.Tensor,
        batch_index: Optional[torch.Tensor] = None,
        give_mean: bool = True,
    ) -> torch.Tensor:
        """
        Provides the tensor of library size from the posterior.

        Parameters
        ----------
        x
            tensor of values with shape ``(batch_size, n_input_genes)``
        y
            tensor of values with shape ``(batch_size, n_input_proteins)``
        batch_index
            tensor of values with shape ``(batch_size, 1)``
        give_mean
            return mean of l or sample from it

        Returns
        -------
        type
            tensor of shape ``(batch_size, 1)``
        """
        if self.log_variational:
            x = torch.log(1 + x)
            y = torch.log(1 + y)
        _, _, ql_m, ql_v, latent, _ = self.encoder(
            torch.cat((x, y), dim=-1), batch_index
        )
        library_gene = latent["l"]
        if give_mean is True:
            return torch.exp(ql_m + 0.5 * ql_v)
        else:
            return library_gene

    def get_sample_dispersion(
        self,
        x: torch.Tensor,
        y: torch.Tensor,
        batch_index: Optional[torch.Tensor] = None,
        label: Optional[torch.Tensor] = None,
        n_samples: int = 1,
    ) -> Tuple[torch.Tensor, torch.Tensor]:
        """
        Returns the tensors of dispersions for genes and proteins.

        Parameters
        ----------
        x
            tensor of values with shape ``(batch_size, n_input_genes)``
        y
            tensor of values with shape ``(batch_size, n_input_proteins)``
        batch_index
            array that indicates which batch the cells belong to with shape ``batch_size``
        label
            tensor of cell-types labels with shape ``(batch_size, n_labels)``
        n_samples
            number of samples

        Returns
        -------
        type
            tensors of dispersions of the negative binomial distribution
        """
        outputs = self.inference(
            x, y, batch_index=batch_index, label=label, n_samples=n_samples
        )
        px_r = outputs["px_"]["r"]
        py_r = outputs["py_"]["r"]
        return px_r, py_r

    def get_reconstruction_loss(
        self,
        x: torch.Tensor,
        y: torch.Tensor,
        px_dict: Dict[str, torch.Tensor],
        py_dict: Dict[str, torch.Tensor],
        pro_batch_mask_minibatch: Optional[torch.Tensor] = None,
    ) -> Tuple[torch.Tensor, torch.Tensor]:
        """Compute reconstruction loss."""
        px_ = px_dict
        py_ = py_dict
        # Reconstruction Loss
        if self.gene_likelihood == "zinb":
            reconst_loss_gene = (
                -ZeroInflatedNegativeBinomial(
                    mu=px_["rate"], theta=px_["r"], zi_logits=px_["dropout"]
                )
                .log_prob(x)
                .sum(dim=-1)
            )
        else:
            reconst_loss_gene = (
                -NegativeBinomial(mu=px_["rate"], theta=px_["r"])
                .log_prob(x)
                .sum(dim=-1)
            )

        py_conditional = NegativeBinomialMixture(
            mu1=py_["rate_back"],
            mu2=py_["rate_fore"],
            theta1=py_["r"],
            mixture_logits=py_["mixing"],
        )
        reconst_loss_protein_full = -py_conditional.log_prob(y)
        if pro_batch_mask_minibatch is not None:
            temp_pro_loss_full = torch.zeros_like(reconst_loss_protein_full)
            temp_pro_loss_full.masked_scatter_(
                pro_batch_mask_minibatch.bool(), reconst_loss_protein_full
            )

            reconst_loss_protein = temp_pro_loss_full.sum(dim=-1)
        else:
            reconst_loss_protein = reconst_loss_protein_full.sum(dim=-1)

        return reconst_loss_gene, reconst_loss_protein

    def inference(
        self,
        x: torch.Tensor,
        y: torch.Tensor,
        batch_index: Optional[torch.Tensor] = None,
        label: Optional[torch.Tensor] = None,
        n_samples=1,
        transform_batch: Optional[int] = None,
    ) -> Dict[str, Union[torch.Tensor, Dict[str, torch.Tensor]]]:
        """
        Internal helper function to compute necessary inference quantities.

        We use the dictionary ``px_`` to contain the parameters of the ZINB/NB for genes.
        The rate refers to the mean of the NB, dropout refers to Bernoulli mixing parameters.
        `scale` refers to the quanity upon which differential expression is performed. For genes,
        this can be viewed as the mean of the underlying gamma distribution.

        We use the dictionary ``py_`` to contain the parameters of the Mixture NB distribution for proteins.
        `rate_fore` refers to foreground mean, while `rate_back` refers to background mean. ``scale`` refers to
        foreground mean adjusted for background probability and scaled to reside in simplex.
        ``back_alpha`` and ``back_beta`` are the posterior parameters for ``rate_back``.  ``fore_scale`` is the scaling
        factor that enforces `rate_fore` > `rate_back`.

        ``px_["r"]`` and ``py_["r"]`` are the inverse dispersion parameters for genes and protein, respectively.
        """
        x_ = x
        y_ = y
        if self.use_observed_lib_size:
            library_gene = x.sum(1).unsqueeze(1)
        if self.log_variational:
            x_ = torch.log(1 + x_)
            y_ = torch.log(1 + y_)

        # Sampling - Encoder gets concatenated genes + proteins
        qz_m, qz_v, ql_m, ql_v, latent, untran_latent = self.encoder(
            torch.cat((x_, y_), dim=-1), batch_index
        )
        z = latent["z"]
        untran_z = untran_latent["z"]
        untran_l = untran_latent["l"]
        if not self.use_observed_lib_size:
            library_gene = latent["l"]

        if n_samples > 1:
            qz_m = qz_m.unsqueeze(0).expand((n_samples, qz_m.size(0), qz_m.size(1)))
            qz_v = qz_v.unsqueeze(0).expand((n_samples, qz_v.size(0), qz_v.size(1)))
            untran_z = Normal(qz_m, qz_v.sqrt()).sample()
            z = self.encoder.z_transformation(untran_z)
            ql_m = ql_m.unsqueeze(0).expand((n_samples, ql_m.size(0), ql_m.size(1)))
            ql_v = ql_v.unsqueeze(0).expand((n_samples, ql_v.size(0), ql_v.size(1)))
            untran_l = Normal(ql_m, ql_v.sqrt()).sample()
            if self.use_observed_lib_size:
                library_gene = library_gene.unsqueeze(0).expand(
                    (n_samples, library_gene.size(0), library_gene.size(1))
                )
            else:
                library_gene = self.encoder.l_transformation(untran_l)

        if self.gene_dispersion == "gene-label":
            # px_r gets transposed - last dimension is nb genes
            px_r = F.linear(one_hot(label, self.n_labels), self.px_r)
        elif self.gene_dispersion == "gene-batch":
            px_r = F.linear(one_hot(batch_index, self.n_batch), self.px_r)
        elif self.gene_dispersion == "gene":
            px_r = self.px_r
        px_r = torch.exp(px_r)

        if self.protein_dispersion == "protein-label":
            # py_r gets transposed - last dimension is n_proteins
            py_r = F.linear(one_hot(label, self.n_labels), self.py_r)
        elif self.protein_dispersion == "protein-batch":
            py_r = F.linear(one_hot(batch_index, self.n_batch), self.py_r)
        elif self.protein_dispersion == "protein":
            py_r = self.py_r
        py_r = torch.exp(py_r)

        # Background regularization
        if self.n_batch > 0:
            py_back_alpha_prior = F.linear(
                one_hot(batch_index, self.n_batch), self.background_pro_alpha
            )
            py_back_beta_prior = F.linear(
                one_hot(batch_index, self.n_batch),
                torch.exp(self.background_pro_log_beta),
            )
        else:
            py_back_alpha_prior = self.background_pro_alpha
            py_back_beta_prior = torch.exp(self.background_pro_log_beta)
        self.back_mean_prior = Normal(py_back_alpha_prior, py_back_beta_prior)

        if transform_batch is not None:
            batch_index = torch.ones_like(batch_index) * transform_batch
        px_, py_, log_pro_back_mean = self.decoder(z, library_gene, batch_index, label)
        px_["r"] = px_r
        py_["r"] = py_r

        return dict(
            px_=px_,
            py_=py_,
            qz_m=qz_m,
            qz_v=qz_v,
            z=z,
            untran_z=untran_z,
            ql_m=ql_m,
            ql_v=ql_v,
            library_gene=library_gene,
            untran_l=untran_l,
            log_pro_back_mean=log_pro_back_mean,
        )

    def forward(
        self,
        x: torch.Tensor,
        y: torch.Tensor,
        local_l_mean_gene: torch.Tensor,
        local_l_var_gene: torch.Tensor,
        batch_index: Optional[torch.Tensor] = None,
        label: Optional[torch.Tensor] = None,
    ) -> Tuple[
        torch.FloatTensor, torch.FloatTensor, torch.FloatTensor, torch.FloatTensor
    ]:
        """
        Returns the reconstruction loss and the Kullback divergences.

        Parameters
        ----------
        x
            tensor of values with shape ``(batch_size, n_input_genes)``
        y
            tensor of values with shape ``(batch_size, n_input_proteins)``
        local_l_mean_gene
            tensor of means of the prior distribution of latent variable l
            with shape ``(batch_size, 1)````
        local_l_var_gene
            tensor of variancess of the prior distribution of latent variable l
            with shape ``(batch_size, 1)``
        batch_index
            array that indicates which batch the cells belong to with shape ``batch_size``
        label
            tensor of cell-types labels with shape (batch_size, n_labels)

        Returns
        -------
        type
            the reconstruction loss and the Kullback divergences
        """
        # Parameters for z latent distribution

        outputs = self.inference(x, y, batch_index, label)
        qz_m = outputs["qz_m"]
        qz_v = outputs["qz_v"]
        ql_m = outputs["ql_m"]
        ql_v = outputs["ql_v"]
        px_ = outputs["px_"]
        py_ = outputs["py_"]

        if self.protein_batch_mask is not None:
            pro_batch_mask_minibatch = torch.zeros_like(y)
            for b in torch.unique(batch_index):
                b_indices = (batch_index == b).reshape(-1)
                pro_batch_mask_minibatch[b_indices] = torch.tensor(
                    self.protein_batch_mask[b.item()].astype(np.float32),
                    device=y.device,
                )
        else:
            pro_batch_mask_minibatch = None

        reconst_loss_gene, reconst_loss_protein = self.get_reconstruction_loss(
            x, y, px_, py_, pro_batch_mask_minibatch
        )

        # KL Divergence
        kl_div_z = kl(Normal(qz_m, torch.sqrt(qz_v)), Normal(0, 1)).sum(dim=1)
        if not self.use_observed_lib_size:
            kl_div_l_gene = kl(
                Normal(ql_m, torch.sqrt(ql_v)),
                Normal(local_l_mean_gene, torch.sqrt(local_l_var_gene)),
            ).sum(dim=1)
        else:
            kl_div_l_gene = 0.0

        kl_div_back_pro_full = kl(
            Normal(py_["back_alpha"], py_["back_beta"]), self.back_mean_prior
        )
        if pro_batch_mask_minibatch is not None:
            kl_div_back_pro = (pro_batch_mask_minibatch * kl_div_back_pro_full).sum(
                dim=1
            )
        else:
            kl_div_back_pro = kl_div_back_pro_full.sum(dim=1)

        return (
            reconst_loss_gene,
            reconst_loss_protein,
            kl_div_z,
            kl_div_l_gene,
            kl_div_back_pro,
        )
