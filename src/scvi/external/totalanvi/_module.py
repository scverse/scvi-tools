"""Main module."""

from collections.abc import Iterable, Sequence
from typing import Literal

import numpy as np
import torch
import torch.nn.functional as F
from torch.distributions import Categorical, Normal
from torch.distributions import kl_divergence as kl
from torch.nn.functional import one_hot

from scvi import REGISTRY_KEYS
from scvi.module import TOTALVAE
from scvi.module._classifier import Classifier
from scvi.module._utils import broadcast_labels
from scvi.module.base import LossOutput, SupervisedModuleClass, auto_move_data
from scvi.nn import Decoder, Encoder


# VAE model
class TOTALANVAE(SupervisedModuleClass, TOTALVAE):
    """Total variational inference for CITE-seq data.

    Implements a combination of scANVI and totalVI model of :cite:p:`GayosoSteier21`.

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
        Dictionary where each key is a batch code, and value is for each protein, whether it was
        observed or not.
    encode_covariates
        Whether to concatenate covariates to expression in encoder
    protein_background_prior_mean
        Array of proteins by batches, the prior initialization for the protein background mean
        (log scale)
    protein_background_prior_scale
        Array of proteins by batches, the prior initialization for the protein background scale
        (log scale)
    use_size_factor_key
        Use size_factor AnnDataField defined by the user as scaling factor in mean of conditional
        distribution. Takes priority over `use_observed_lib_size`.
    use_observed_lib_size
        Use observed library size for RNA as scaling factor in mean of conditional distribution
    extra_payload_autotune
        If ``True``, will return extra matrices in the loss output to be used during autotune
    library_log_means
        1 x n_batch array of means of the log library sizes. Parameterizes prior on library size if
        not using observed library size.
    library_log_vars
        1 x n_batch array of variances of the log library sizes. Parameterizes prior on library
        size if not using observed library size.
    use_batch_norm
        Whether to use batch norm in layers.
    use_layer_norm
        Whether to use layer norm in layers.
    extra_encoder_kwargs
        Extra keyword arguments passed into :class:`~scvi.nn.EncoderTOTALVI`.
    extra_decoder_kwargs
        Extra keyword arguments passed into :class:`~scvi.nn.DecoderTOTALVI`.
    """

    def __init__(
        self,
        n_input_genes: int,
        n_input_proteins: int,
        n_batch: int = 1,
        n_labels: int = 1,
        n_hidden: int = 256,
        n_latent: int = 20,
        n_layers_encoder: int = 2,
        n_layers_decoder: int = 1,
        n_continuous_cov: int = 0,
        n_cats_per_cov: Iterable[int] | None = None,
        dropout_rate_decoder: float = 0.2,
        dropout_rate_encoder: float = 0.2,
        gene_dispersion: Literal["gene", "gene-batch", "gene-label"] = "gene",
        protein_dispersion: Literal["protein", "protein-batch", "protein-label"] = "protein",
        log_variational: bool = True,
        gene_likelihood: Literal["zinb", "nb"] = "nb",
        latent_distribution: Literal["normal", "ln"] = "normal",
        protein_batch_mask: dict[str | int, np.ndarray] = None,
        encode_covariates: bool = True,
        protein_background_prior_mean: np.ndarray | None = None,
        protein_background_prior_scale: np.ndarray | None = None,
        use_size_factor_key: bool = False,
        use_observed_lib_size: bool = True,
        extra_payload_autotune: bool = False,
        library_log_means: np.ndarray | None = None,
        library_log_vars: np.ndarray | None = None,
        n_panel: int | None = None,
        panel_key: str = REGISTRY_KEYS.BATCH_KEY,
        use_batch_norm: Literal["encoder", "decoder", "none", "both"] = "both",
        use_layer_norm: Literal["encoder", "decoder", "none", "both"] = "none",
        extra_encoder_kwargs: dict | None = None,
        extra_decoder_kwargs: dict | None = None,
        y_prior=None,
        labels_groups: Sequence[int] = None,
        use_labels_groups: bool = False,
        linear_classifier: bool = False,
        classifier_parameters: dict | None = None,
    ):
        super().__init__(
            n_input_genes,
            n_input_proteins,
            n_batch=n_batch,
            n_panel=n_panel,
            panel_key=panel_key,
            n_labels=n_labels,
            n_hidden=n_hidden,
            n_latent=n_latent,
            n_layers_encoder=n_layers_encoder,
            n_layers_decoder=n_layers_decoder,
            n_continuous_cov=n_continuous_cov,
            n_cats_per_cov=n_cats_per_cov,
            dropout_rate_decoder=dropout_rate_decoder,
            dropout_rate_encoder=dropout_rate_encoder,
            gene_dispersion=gene_dispersion,
            protein_dispersion=protein_dispersion,
            log_variational=log_variational,
            gene_likelihood=gene_likelihood,
            latent_distribution=latent_distribution,
            protein_batch_mask=protein_batch_mask,
            encode_covariates=encode_covariates,
            protein_background_prior_mean=protein_background_prior_mean,
            protein_background_prior_scale=protein_background_prior_scale,
            use_size_factor_key=use_size_factor_key,
            use_observed_lib_size=use_observed_lib_size,
            extra_payload_autotune=extra_payload_autotune,
            library_log_means=library_log_means,
            library_log_vars=library_log_vars,
            use_batch_norm=use_batch_norm,
            use_layer_norm=use_layer_norm,
            extra_encoder_kwargs=extra_encoder_kwargs,
            extra_decoder_kwargs=extra_decoder_kwargs,
        )

        self.n_labels = n_labels
        classifier_parameters = classifier_parameters or {}
        use_batch_norm_encoder = use_batch_norm == "encoder" or use_batch_norm == "both"
        use_batch_norm_decoder = use_batch_norm == "decoder" or use_batch_norm == "both"
        use_layer_norm_encoder = use_layer_norm == "encoder" or use_layer_norm == "both"
        use_layer_norm_decoder = use_layer_norm == "decoder" or use_layer_norm == "both"

        # z encoder goes from the n_input-dimensional data to an n_latent-d
        # latent space representation

        self.n_labels = n_labels
        cls_parameters = {
            "n_layers": 0 if linear_classifier else n_layers_encoder,
            "n_hidden": 0 if linear_classifier else n_hidden,
            "dropout_rate": dropout_rate_encoder,
            "logits": True,
        }
        cls_parameters.update(classifier_parameters)
        self.classifier = Classifier(
            n_latent,
            n_labels=n_labels,
            use_batch_norm=use_batch_norm_encoder,
            use_layer_norm=use_layer_norm_encoder,
            **cls_parameters,
        )

        self.encoder_z2_z1 = Encoder(
            n_latent,
            n_latent,
            n_cat_list=[self.n_labels],
            n_layers=n_layers_encoder,
            n_hidden=n_hidden,
            dropout_rate=0.05,
            use_batch_norm=use_batch_norm_encoder,
            use_layer_norm=use_layer_norm_encoder,
            return_dist=True,
        )

        self.decoder_z1_z2 = Decoder(
            n_latent,
            n_latent,
            n_cat_list=[self.n_labels],
            n_layers=n_layers_decoder,
            n_hidden=n_hidden,
            use_batch_norm=use_batch_norm_decoder,
            use_layer_norm=use_layer_norm_decoder,
        )

        self.y_prior = torch.nn.Parameter(
            y_prior if y_prior is not None else (1 / n_labels) * torch.ones(1, n_labels),
            requires_grad=False,
        )
        self.use_labels_groups = use_labels_groups
        self.labels_groups = np.array(labels_groups) if labels_groups is not None else None
        if self.use_labels_groups:
            if labels_groups is None:
                raise ValueError("Specify label groups")
            unique_groups = np.unique(self.labels_groups)
            self.n_groups = len(unique_groups)
            if not (unique_groups == np.arange(self.n_groups)).all():
                raise ValueError()
            self.classifier_groups = Classifier(
                n_latent, n_hidden, self.n_groups, n_layers_encoder, dropout_rate_encoder
            )
            self.groups_index = torch.nn.ParameterList(
                [
                    torch.nn.Parameter(
                        torch.tensor(
                            (self.labels_groups == i).astype(np.uint8),
                            dtype=torch.uint8,
                        ),
                        requires_grad=False,
                    )
                    for i in range(self.n_groups)
                ]
            )

    @auto_move_data
    def classify(
        self,
        x: torch.Tensor,
        y: torch.Tensor,
        batch_index: torch.Tensor | None = None,
        cont_covs: torch.Tensor | None = None,
        cat_covs: torch.Tensor | None = None,
        use_posterior_mean: bool = True,
    ) -> torch.Tensor:
        """Forward pass through the encoder and classifier.

        Parameters
        ----------
        x
            Tensor of shape ``(n_obs, n_genes)``.
        y
            Tensor of shape ``(n_obs, n_proteins)``.
        batch_index
            Tensor of shape ``(n_obs,)`` denoting batch indices.
        cont_covs
            Tensor of shape ``(n_obs, n_continuous_covariates)``.
        cat_covs
            Tensor of shape ``(n_obs, n_categorical_covariates)``.
        use_posterior_mean
            Whether to use the posterior mean of the latent distribution for
            classification.

        Returns
        -------
        Tensor of shape ``(n_obs, n_labels)`` denoting logit scores per label.
        Before v1.1, this method by default returned probabilities per label,
        see #2301 for more details.
        """
        x_ = x / (1 + x.mean(1, keepdim=True))
        y_ = y / (1 + y.mean(1, keepdim=True))

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
            categorical_input = ()
        qz, _, _, _ = self.encoder(encoder_input, batch_index, *categorical_input)

        z = qz.loc if use_posterior_mean else qz.rsample()

        return self.classify_helper(z)

    def loss(
        self,
        tensors,
        inference_outputs,
        generative_outputs,
        pro_recons_weight=1.0,  # double check these defaults
        kl_weight=1.0,
        labelled_tensors=None,
        classification_ratio=None,
    ) -> tuple[torch.FloatTensor, torch.FloatTensor, torch.FloatTensor, torch.FloatTensor]:
        """Returns the reconstruction loss and the Kullback divergences.

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
        qz1 = inference_outputs["qz"]
        z1 = inference_outputs["z"]
        ql = inference_outputs["ql"]
        px_ = generative_outputs["px_"]
        py_ = generative_outputs["py_"]
        per_batch_efficiency = generative_outputs["per_batch_efficiency"]

        x = tensors[REGISTRY_KEYS.X_KEY]
        batch_index = tensors[REGISTRY_KEYS.BATCH_KEY]
        panel_index = tensors[self.panel_key]
        label_index = tensors[REGISTRY_KEYS.LABELS_KEY]
        y = tensors[REGISTRY_KEYS.PROTEIN_EXP_KEY]

        # Enumerate choices of label
        labels_, z1s = broadcast_labels(z1, n_broadcast=self.n_labels)
        qz2, z2 = self.encoder_z2_z1(z1s, labels_)
        pz1_m, pz1_v = self.decoder_z1_z2(z2, labels_)

        # KL Divergence
        mean = torch.zeros_like(qz2.loc)
        scale = torch.ones_like(qz2.scale)

        kl_divergence_z2 = kl(qz2, Normal(mean, scale)).sum(dim=1)
        loss_z1_unweight = -Normal(pz1_m, torch.sqrt(pz1_v)).log_prob(z1s).sum(dim=-1)
        loss_z1_weight = qz1.log_prob(z1).sum(dim=-1)

        if not self.use_observed_lib_size:
            n_batch = self.library_log_means.shape[1]
            local_library_log_means = F.linear(
                one_hot(batch_index.squeeze(-1), n_batch).float(), self.library_log_means
            )
            local_library_log_vars = F.linear(
                one_hot(batch_index.squeeze(-1), n_batch).float(), self.library_log_vars
            )
            kl_div_l_gene = kl(
                ql,
                Normal(local_library_log_means, torch.sqrt(local_library_log_vars)),
            ).sum(dim=1)
        else:
            kl_div_l_gene = 0.0

        if self.protein_batch_mask is not None:
            pro_batch_mask_minibatch = torch.zeros_like(y)
            for b in torch.unique(panel_index):
                b_indices = (panel_index == b).reshape(-1)
                pro_batch_mask_minibatch[b_indices] = torch.tensor(
                    self.protein_batch_mask[str(int(b.item()))].astype(np.float32),
                    device=y.device,
                )
        else:
            pro_batch_mask_minibatch = None

        reconst_loss_gene, reconst_loss_protein = self.get_reconstruction_loss(
            x, y, px_, py_, pro_batch_mask_minibatch, per_batch_efficiency
        )

        kl_div_back_pro_full = kl(
            Normal(py_["back_alpha"], py_["back_beta"]), inference_outputs["back_mean_prior"]
        )
        lkl_back_pro_full = -torch.distributions.LogNormal(
            torch.tensor([0.0]).to(x.device), torch.tensor([1.0]).to(x.device)
        ).log_prob(per_batch_efficiency)
        lkl_protein_expressed = -1e-3 * torch.distributions.Bernoulli(
            logits=py_["mixing"]
        ).log_prob(torch.ones_like(py_["mixing"]))
        if pro_batch_mask_minibatch is not None:
            kl_div_back_pro = pro_batch_mask_minibatch.bool() * kl_div_back_pro_full
            kl_div_back_pro = (
                kl_div_back_pro.sum(dim=1)
                + lkl_back_pro_full.sum(dim=1)
                + lkl_protein_expressed.sum(dim=1)
            )
        else:
            kl_div_back_pro = (
                kl_div_back_pro_full.sum(dim=1)
                + lkl_back_pro_full.sum(dim=1)
                + lkl_protein_expressed.sum(dim=1)
            )

        reconst_loss = reconst_loss_gene + kl_weight * pro_recons_weight * reconst_loss_protein
        reconst_losses = {
            "reconst_loss_gene": reconst_loss_gene,
            "reconst_loss_protein": reconst_loss_protein,
        }
        kl_locals = {
            # "kl_div_l_gene": kl_div_l_gene,
            "kl_div_back_pro": kl_div_back_pro,
        }

        probs = self.classifier(z1)
        if self.classifier.logits:
            probs = F.softmax(probs, dim=-1)

        reconst_loss += loss_z1_weight + (
            (loss_z1_unweight).view(self.n_labels, -1).t() * probs
        ).sum(dim=1)

        torch.stack(
            [
                torch.nn.functional.one_hot(y_i, self.n_labels)
                if y_i < self.n_labels
                else probs[i, :]
                for i, y_i in enumerate(label_index.squeeze(-1))
            ]
        )
        kl_divergence = (kl_divergence_z2.view(self.n_labels, -1).t() * probs).sum(dim=1)
        kl_locals["kl_divergence"] = kl_divergence
        kl_divergence_class = kl(
            Categorical(probs=probs),
            Categorical(probs=self.y_prior.repeat(probs.size(0), 1)),
        )
        kl_locals["kl_divergence_class"] = kl_divergence_class

        loss = torch.mean(
            reconst_loss
            + kl_weight * kl_divergence
            + kl_div_l_gene
            + kl_weight * kl_div_back_pro
            + kl_weight * kl_divergence_class
        )

        # a payload to be used during autotune
        if self.extra_payload_autotune:
            extra_metrics_payload = {
                "z": inference_outputs["z"],
                "batch": tensors[REGISTRY_KEYS.BATCH_KEY],
                "labels": tensors[REGISTRY_KEYS.LABELS_KEY],
            }
        else:
            extra_metrics_payload = {}

        if labelled_tensors is not None:
            ce_loss, true_labels, logits = self.classification_loss(labelled_tensors)
            loss += ce_loss * classification_ratio
            return LossOutput(
                loss=loss,
                reconstruction_loss=reconst_losses,
                kl_local=kl_locals,
                classification_loss=ce_loss,
                true_labels=true_labels,
                logits=logits,
                extra_metrics=extra_metrics_payload,
            )
        return LossOutput(
            loss=loss,
            reconstruction_loss=reconst_losses,
            kl_local=kl_locals,
            extra_metrics=extra_metrics_payload,
        )
