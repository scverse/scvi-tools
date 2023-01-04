"""Main module."""
from typing import Dict, Iterable, Literal, Optional, Tuple, Union

import numpy as np
import torch
import torch.nn.functional as F
from torch.distributions import Normal
from torch.distributions import kl_divergence as kl

from scvi import REGISTRY_KEYS
from scvi.autotune._types import Tunable
from scvi.distributions import (
    NegativeBinomial,
    NegativeBinomialMixture,
    ZeroInflatedNegativeBinomial,
)
from scvi.module.base import BaseModuleClass, LossOutput, auto_move_data
from scvi.nn import DecoderTOTALVI, EncoderTOTALVI, one_hot

torch.backends.cudnn.benchmark = True


# VAE model
class TOTALVAE(BaseModuleClass):
    """
    Total variational inference for CITE-seq data.

    Implements the totalVI model of :cite:p:`GayosoSteier21`.

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
    n_continuous_cov
        Number of continuous covarites
    n_cats_per_cov
        Number of categories for each extra categorical covariate
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
    use_size_factor_key
        Use size_factor AnnDataField defined by the user as scaling factor in mean of conditional distribution.
        Takes priority over `use_observed_lib_size`.
    use_observed_lib_size
        Use observed library size for RNA as scaling factor in mean of conditional distribution
    library_log_means
        1 x n_batch array of means of the log library sizes. Parameterizes prior on library size if
        not using observed library size.
    library_log_vars
        1 x n_batch array of variances of the log library sizes. Parameterizes prior on library size if
        not using observed library size.
    """

    def __init__(
        self,
        n_input_genes: int,
        n_input_proteins: int,
        n_batch: int = 0,
        n_labels: int = 0,
        n_hidden: Tunable[int] = 256,
        n_latent: Tunable[int] = 20,
        n_layers_encoder: Tunable[int] = 2,
        n_layers_decoder: Tunable[int] = 1,
        n_continuous_cov: int = 0,
        n_cats_per_cov: Optional[Iterable[int]] = None,
        dropout_rate_decoder: Tunable[float] = 0.2,
        dropout_rate_encoder: Tunable[float] = 0.2,
        gene_dispersion: Tunable[Literal["gene", "gene-batch", "gene-label"]] = "gene",
        protein_dispersion: Tunable[
            Literal["protein", "protein-batch", "protein-label"]
        ] = "protein",
        log_variational: bool = True,
        gene_likelihood: Tunable[Literal["zinb", "nb"]] = "nb",
        latent_distribution: Tunable[Literal["normal", "ln"]] = "normal",
        protein_batch_mask: Dict[Union[str, int], np.ndarray] = None,
        encode_covariates: bool = True,
        protein_background_prior_mean: Optional[np.ndarray] = None,
        protein_background_prior_scale: Optional[np.ndarray] = None,
        use_size_factor_key: bool = False,
        use_observed_lib_size: bool = True,
        library_log_means: Optional[np.ndarray] = None,
        library_log_vars: Optional[np.ndarray] = None,
        use_batch_norm: Tunable[Literal["encoder", "decoder", "none", "both"]] = "both",
        use_layer_norm: Tunable[Literal["encoder", "decoder", "none", "both"]] = "none",
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
        self.encode_covariates = encode_covariates
        self.use_size_factor_key = use_size_factor_key
        self.use_observed_lib_size = use_size_factor_key or use_observed_lib_size
        if not self.use_observed_lib_size:
            if library_log_means is None or library_log_means is None:
                raise ValueError(
                    "If not using observed_lib_size, "
                    "must provide library_log_means and library_log_vars."
                )

            self.register_buffer(
                "library_log_means", torch.from_numpy(library_log_means).float()
            )
            self.register_buffer(
                "library_log_vars", torch.from_numpy(library_log_vars).float()
            )

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
        n_input = n_input_genes + self.n_input_proteins
        n_input_encoder = n_input + n_continuous_cov * encode_covariates
        cat_list = [n_batch] + list([] if n_cats_per_cov is None else n_cats_per_cov)
        encoder_cat_list = cat_list if encode_covariates else None
        self.encoder = EncoderTOTALVI(
            n_input_encoder,
            n_latent,
            n_layers=n_layers_encoder,
            n_cat_list=encoder_cat_list,
            n_hidden=n_hidden,
            dropout_rate=dropout_rate_encoder,
            distribution=latent_distribution,
            use_batch_norm=use_batch_norm_encoder,
            use_layer_norm=use_layer_norm_encoder,
        )
        self.decoder = DecoderTOTALVI(
            n_latent + n_continuous_cov,
            n_input_genes,
            self.n_input_proteins,
            n_layers=n_layers_decoder,
            n_cat_list=cat_list,
            n_hidden=n_hidden,
            dropout_rate=dropout_rate_decoder,
            use_batch_norm=use_batch_norm_decoder,
            use_layer_norm=use_layer_norm_decoder,
            scale_activation="softplus" if use_size_factor_key else "softmax",
        )

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

    def _get_inference_input(self, tensors):
        x = tensors[REGISTRY_KEYS.X_KEY]
        y = tensors[REGISTRY_KEYS.PROTEIN_EXP_KEY]
        batch_index = tensors[REGISTRY_KEYS.BATCH_KEY]

        cont_key = REGISTRY_KEYS.CONT_COVS_KEY
        cont_covs = tensors[cont_key] if cont_key in tensors.keys() else None

        cat_key = REGISTRY_KEYS.CAT_COVS_KEY
        cat_covs = tensors[cat_key] if cat_key in tensors.keys() else None

        input_dict = dict(
            x=x, y=y, batch_index=batch_index, cat_covs=cat_covs, cont_covs=cont_covs
        )
        return input_dict

    def _get_generative_input(self, tensors, inference_outputs):
        z = inference_outputs["z"]
        library_gene = inference_outputs["library_gene"]
        batch_index = tensors[REGISTRY_KEYS.BATCH_KEY]
        label = tensors[REGISTRY_KEYS.LABELS_KEY]

        cont_key = REGISTRY_KEYS.CONT_COVS_KEY
        cont_covs = tensors[cont_key] if cont_key in tensors.keys() else None

        cat_key = REGISTRY_KEYS.CAT_COVS_KEY
        cat_covs = tensors[cat_key] if cat_key in tensors.keys() else None

        size_factor_key = REGISTRY_KEYS.SIZE_FACTOR_KEY
        size_factor = (
            tensors[size_factor_key] if size_factor_key in tensors.keys() else None
        )

        return dict(
            z=z,
            library_gene=library_gene,
            batch_index=batch_index,
            label=label,
            cat_covs=cat_covs,
            cont_covs=cont_covs,
            size_factor=size_factor,
        )

    @auto_move_data
    def generative(
        self,
        z: torch.Tensor,
        library_gene: torch.Tensor,
        batch_index: torch.Tensor,
        label: torch.Tensor,
        cont_covs=None,
        cat_covs=None,
        size_factor=None,
        transform_batch: Optional[int] = None,
    ) -> Dict[str, Union[torch.Tensor, Dict[str, torch.Tensor]]]:
        """Run the generative step."""
        if cont_covs is None:
            decoder_input = z
        elif z.dim() != cont_covs.dim():
            decoder_input = torch.cat(
                [z, cont_covs.unsqueeze(0).expand(z.size(0), -1, -1)], dim=-1
            )
        else:
            decoder_input = torch.cat([z, cont_covs], dim=-1)

        if cat_covs is not None:
            categorical_input = torch.split(cat_covs, 1, dim=1)
        else:
            categorical_input = tuple()

        if transform_batch is not None:
            batch_index = torch.ones_like(batch_index) * transform_batch

        if not self.use_size_factor_key:
            size_factor = library_gene

        px_, py_, log_pro_back_mean = self.decoder(
            decoder_input, size_factor, batch_index, *categorical_input
        )

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

        px_["r"] = px_r
        py_["r"] = py_r
        return dict(
            px_=px_,
            py_=py_,
            log_pro_back_mean=log_pro_back_mean,
        )

    @auto_move_data
    def inference(
        self,
        x: torch.Tensor,
        y: torch.Tensor,
        batch_index: Optional[torch.Tensor] = None,
        label: Optional[torch.Tensor] = None,
        n_samples=1,
        cont_covs=None,
        cat_covs=None,
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

        Parameters
        ----------
        x
            tensor of values with shape ``(batch_size, n_input_genes)``
        y
            tensor of values with shape ``(batch_size, n_input_proteins)``
        batch_index
            array that indicates which batch the cells belong to with shape ``batch_size``
        label
            tensor of cell-types labels with shape (batch_size, n_labels)
        n_samples
            Number of samples to sample from approximate posterior
        cont_covs
            Continuous covariates to condition on
        cat_covs
            Categorical covariates to condition on
        """
        x_ = x
        y_ = y
        if self.use_observed_lib_size:
            library_gene = x.sum(1).unsqueeze(1)
        if self.log_variational:
            x_ = torch.log(1 + x_)
            y_ = torch.log(1 + y_)

        if cont_covs is not None and self.encode_covariates is True:
            encoder_input = torch.cat((x_, y_, cont_covs), dim=-1)
        else:
            encoder_input = torch.cat((x_, y_), dim=-1)
        if cat_covs is not None and self.encode_covariates is True:
            categorical_input = torch.split(cat_covs, 1, dim=1)
        else:
            categorical_input = tuple()
        qz, ql, latent, untran_latent = self.encoder(
            encoder_input, batch_index, *categorical_input
        )

        z = latent["z"]
        untran_z = untran_latent["z"]
        untran_l = untran_latent["l"]
        if not self.use_observed_lib_size:
            library_gene = latent["l"]

        if n_samples > 1:
            untran_z = qz.sample((n_samples,))
            z = self.encoder.z_transformation(untran_z)

            untran_l = ql.sample((n_samples,))
            if self.use_observed_lib_size:
                library_gene = library_gene.unsqueeze(0).expand(
                    (n_samples, library_gene.size(0), library_gene.size(1))
                )
            else:
                library_gene = self.encoder.l_transformation(untran_l)

        # Background regularization
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

        return dict(
            qz=qz,
            z=z,
            untran_z=untran_z,
            ql=ql,
            library_gene=library_gene,
            untran_l=untran_l,
        )

    def loss(
        self,
        tensors,
        inference_outputs,
        generative_outputs,
        pro_recons_weight=1.0,  # double check these defaults
        kl_weight=1.0,
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
        batch_index
            array that indicates which batch the cells belong to with shape ``batch_size``
        label
            tensor of cell-types labels with shape (batch_size, n_labels)

        Returns
        -------
        type
            the reconstruction loss and the Kullback divergences
        """
        qz = inference_outputs["qz"]
        ql = inference_outputs["ql"]
        px_ = generative_outputs["px_"]
        py_ = generative_outputs["py_"]

        x = tensors[REGISTRY_KEYS.X_KEY]
        batch_index = tensors[REGISTRY_KEYS.BATCH_KEY]
        y = tensors[REGISTRY_KEYS.PROTEIN_EXP_KEY]

        if self.protein_batch_mask is not None:
            pro_batch_mask_minibatch = torch.zeros_like(y)
            for b in torch.unique(batch_index):
                b_indices = (batch_index == b).reshape(-1)
                pro_batch_mask_minibatch[b_indices] = torch.tensor(
                    self.protein_batch_mask[str(int(b.item()))].astype(np.float32),
                    device=y.device,
                )
        else:
            pro_batch_mask_minibatch = None

        reconst_loss_gene, reconst_loss_protein = self.get_reconstruction_loss(
            x, y, px_, py_, pro_batch_mask_minibatch
        )

        # KL Divergence
        kl_div_z = kl(qz, Normal(0, 1)).sum(dim=1)
        if not self.use_observed_lib_size:
            n_batch = self.library_log_means.shape[1]
            local_library_log_means = F.linear(
                one_hot(batch_index, n_batch), self.library_log_means
            )
            local_library_log_vars = F.linear(
                one_hot(batch_index, n_batch), self.library_log_vars
            )
            kl_div_l_gene = kl(
                ql,
                Normal(local_library_log_means, torch.sqrt(local_library_log_vars)),
            ).sum(dim=1)
        else:
            kl_div_l_gene = 0.0

        kl_div_back_pro_full = kl(
            Normal(py_["back_alpha"], py_["back_beta"]), self.back_mean_prior
        )
        if pro_batch_mask_minibatch is not None:
            kl_div_back_pro = torch.zeros_like(kl_div_back_pro_full)
            kl_div_back_pro.masked_scatter_(
                pro_batch_mask_minibatch.bool(), kl_div_back_pro_full
            )
            kl_div_back_pro = kl_div_back_pro.sum(dim=1)
        else:
            kl_div_back_pro = kl_div_back_pro_full.sum(dim=1)
        loss = torch.mean(
            reconst_loss_gene
            + pro_recons_weight * reconst_loss_protein
            + kl_weight * kl_div_z
            + kl_div_l_gene
            + kl_weight * kl_div_back_pro
        )

        reconst_losses = dict(
            reconst_loss_gene=reconst_loss_gene,
            reconst_loss_protein=reconst_loss_protein,
        )
        kl_local = dict(
            kl_div_z=kl_div_z,
            kl_div_l_gene=kl_div_l_gene,
            kl_div_back_pro=kl_div_back_pro,
        )

        return LossOutput(
            loss=loss, reconstruction_loss=reconst_losses, kl_local=kl_local
        )

    @torch.inference_mode()
    def sample(self, tensors, n_samples=1):
        """Sample from the generative model."""
        inference_kwargs = dict(n_samples=n_samples)
        with torch.inference_mode():
            inference_outputs, generative_outputs, = self.forward(
                tensors,
                inference_kwargs=inference_kwargs,
                compute_loss=False,
            )

        px_ = generative_outputs["px_"]
        py_ = generative_outputs["py_"]

        rna_dist = NegativeBinomial(mu=px_["rate"], theta=px_["r"])
        protein_dist = NegativeBinomialMixture(
            mu1=py_["rate_back"],
            mu2=py_["rate_fore"],
            theta1=py_["r"],
            mixture_logits=py_["mixing"],
        )
        rna_sample = rna_dist.sample().cpu()
        protein_sample = protein_dist.sample().cpu()

        return rna_sample, protein_sample

    @torch.inference_mode()
    @auto_move_data
    def marginal_ll(self, tensors, n_mc_samples):
        """Computes the marginal log likelihood of the data under the model."""
        x = tensors[REGISTRY_KEYS.X_KEY]
        batch_index = tensors[REGISTRY_KEYS.BATCH_KEY]
        to_sum = torch.zeros(x.size()[0], n_mc_samples)

        for i in range(n_mc_samples):
            # Distribution parameters and sampled variables
            inference_outputs, generative_outputs, losses = self.forward(tensors)
            # outputs = self.module.inference(x, y, batch_index, labels)
            qz = inference_outputs["qz"]
            ql = inference_outputs["ql"]
            py_ = generative_outputs["py_"]
            log_library = inference_outputs["untran_l"]
            # really need not softmax transformed random variable
            z = inference_outputs["untran_z"]
            log_pro_back_mean = generative_outputs["log_pro_back_mean"]

            # Reconstruction Loss
            reconst_loss = losses.reconstruction_loss
            reconst_loss_gene = reconst_loss["reconst_loss_gene"]
            reconst_loss_protein = reconst_loss["reconst_loss_protein"]

            # Log-probabilities
            log_prob_sum = torch.zeros(qz.loc.shape[0]).to(self.device)

            if not self.use_observed_lib_size:
                n_batch = self.library_log_means.shape[1]
                local_library_log_means = F.linear(
                    one_hot(batch_index, n_batch), self.library_log_means
                )
                local_library_log_vars = F.linear(
                    one_hot(batch_index, n_batch), self.library_log_vars
                )
                p_l_gene = (
                    Normal(local_library_log_means, local_library_log_vars.sqrt())
                    .log_prob(log_library)
                    .sum(dim=-1)
                )
                q_l_x = ql.log_prob(log_library).sum(dim=-1)

                log_prob_sum += p_l_gene - q_l_x

            p_z = Normal(0, 1).log_prob(z).sum(dim=-1)
            p_mu_back = self.back_mean_prior.log_prob(log_pro_back_mean).sum(dim=-1)
            p_xy_zl = -(reconst_loss_gene + reconst_loss_protein)
            q_z_x = qz.log_prob(z).sum(dim=-1)
            q_mu_back = (
                Normal(py_["back_alpha"], py_["back_beta"])
                .log_prob(log_pro_back_mean)
                .sum(dim=-1)
            )
            log_prob_sum += p_z + p_mu_back + p_xy_zl - q_z_x - q_mu_back

            to_sum[:, i] = log_prob_sum

        batch_log_lkl = torch.logsumexp(to_sum, dim=-1) - np.log(n_mc_samples)
        log_lkl = torch.sum(batch_log_lkl).item()
        return log_lkl
