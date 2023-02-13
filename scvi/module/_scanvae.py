from typing import Iterable, Literal, Optional, Sequence

import numpy as np
import torch
from torch.distributions import Categorical, Normal
from torch.distributions import kl_divergence as kl
from torch.nn import functional as F

from scvi import REGISTRY_KEYS
from scvi.autotune._types import Tunable
from scvi.module.base import LossOutput, auto_move_data
from scvi.nn import Decoder, Encoder

from ._classifier import Classifier
from ._utils import broadcast_labels
from ._vae import VAE


class SCANVAE(VAE):
    """
    Single-cell annotation using variational inference.

    This is an implementation of the scANVI model described in :cite:p:`Xu21`,
    inspired from M1 + M2 model, as described in (https://arxiv.org/pdf/1406.5298.pdf).

    Parameters
    ----------
    n_input
        Number of input genes
    n_batch
        Number of batches
    n_labels
        Number of labels
    n_hidden
        Number of nodes per hidden layer
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
    dispersion
        One of the following

        * ``'gene'`` - dispersion parameter of NB is constant per gene across cells
        * ``'gene-batch'`` - dispersion can differ between different batches
        * ``'gene-label'`` - dispersion can differ between different labels
        * ``'gene-cell'`` - dispersion can differ for every gene in every cell
    log_variational
        Log(data+1) prior to encoding for numerical stability. Not normalization.
    gene_likelihood
        One of:

        * ``'nb'`` - Negative binomial distribution
        * ``'zinb'`` - Zero-inflated negative binomial distribution
    y_prior
        If None, initialized to uniform probability over cell types
    labels_groups
        Label group designations
    use_labels_groups
        Whether to use the label groups
    use_batch_norm
        Whether to use batch norm in layers
    use_layer_norm
        Whether to use layer norm in layers
    **vae_kwargs
        Keyword args for :class:`~scvi.module.VAE`
    """

    def __init__(
        self,
        n_input: int,
        n_batch: int = 0,
        n_labels: int = 0,
        n_hidden: Tunable[int] = 128,
        n_latent: Tunable[int] = 10,
        n_layers: Tunable[int] = 1,
        n_continuous_cov: int = 0,
        n_cats_per_cov: Optional[Iterable[int]] = None,
        dropout_rate: Tunable[float] = 0.1,
        dispersion: Tunable[
            Literal["gene", "gene-batch", "gene-label", "gene-cell"]
        ] = "gene",
        log_variational: Tunable[bool] = True,
        gene_likelihood: Tunable[Literal["zinb", "nb"]] = "zinb",
        y_prior=None,
        labels_groups: Sequence[int] = None,
        use_labels_groups: bool = False,
        classifier_parameters: Optional[dict] = None,
        use_batch_norm: Tunable[Literal["encoder", "decoder", "none", "both"]] = "both",
        use_layer_norm: Tunable[Literal["encoder", "decoder", "none", "both"]] = "none",
        **vae_kwargs,
    ):
        super().__init__(
            n_input,
            n_hidden=n_hidden,
            n_latent=n_latent,
            n_layers=n_layers,
            n_continuous_cov=n_continuous_cov,
            n_cats_per_cov=n_cats_per_cov,
            dropout_rate=dropout_rate,
            n_batch=n_batch,
            dispersion=dispersion,
            log_variational=log_variational,
            gene_likelihood=gene_likelihood,
            use_batch_norm=use_batch_norm,
            use_layer_norm=use_layer_norm,
            **vae_kwargs,
        )

        classifier_parameters = classifier_parameters or dict()
        use_batch_norm_encoder = use_batch_norm == "encoder" or use_batch_norm == "both"
        use_batch_norm_decoder = use_batch_norm == "decoder" or use_batch_norm == "both"
        use_layer_norm_encoder = use_layer_norm == "encoder" or use_layer_norm == "both"
        use_layer_norm_decoder = use_layer_norm == "decoder" or use_layer_norm == "both"

        self.n_labels = n_labels
        # Classifier takes n_latent as input
        cls_parameters = {
            "n_layers": n_layers,
            "n_hidden": n_hidden,
            "dropout_rate": dropout_rate,
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
            n_layers=n_layers,
            n_hidden=n_hidden,
            dropout_rate=dropout_rate,
            use_batch_norm=use_batch_norm_encoder,
            use_layer_norm=use_layer_norm_encoder,
            return_dist=True,
        )

        self.decoder_z1_z2 = Decoder(
            n_latent,
            n_latent,
            n_cat_list=[self.n_labels],
            n_layers=n_layers,
            n_hidden=n_hidden,
            use_batch_norm=use_batch_norm_decoder,
            use_layer_norm=use_layer_norm_decoder,
        )

        self.y_prior = torch.nn.Parameter(
            y_prior
            if y_prior is not None
            else (1 / n_labels) * torch.ones(1, n_labels),
            requires_grad=False,
        )
        self.use_labels_groups = use_labels_groups
        self.labels_groups = (
            np.array(labels_groups) if labels_groups is not None else None
        )
        if self.use_labels_groups:
            if labels_groups is None:
                raise ValueError("Specify label groups")
            unique_groups = np.unique(self.labels_groups)
            self.n_groups = len(unique_groups)
            if not (unique_groups == np.arange(self.n_groups)).all():
                raise ValueError()
            self.classifier_groups = Classifier(
                n_latent, n_hidden, self.n_groups, n_layers, dropout_rate
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
    def classify(self, x, batch_index=None, cont_covs=None, cat_covs=None):
        """Classify cells into cell types."""
        if self.log_variational:
            x = torch.log(1 + x)

        if cont_covs is not None and self.encode_covariates:
            encoder_input = torch.cat((x, cont_covs), dim=-1)
        else:
            encoder_input = x
        if cat_covs is not None and self.encode_covariates:
            categorical_input = torch.split(cat_covs, 1, dim=1)
        else:
            categorical_input = tuple()

        qz, z = self.z_encoder(encoder_input, batch_index, *categorical_input)
        # We classify using the inferred mean parameter of z_1 in the latent space
        z = qz.loc
        if self.use_labels_groups:
            w_g = self.classifier_groups(z)
            unw_y = self.classifier(z)
            w_y = torch.zeros_like(unw_y)
            for i, group_index in enumerate(self.groups_index):
                unw_y_g = unw_y[:, group_index]
                w_y[:, group_index] = unw_y_g / (
                    unw_y_g.sum(dim=-1, keepdim=True) + 1e-8
                )
                w_y[:, group_index] *= w_g[:, [i]]
        else:
            w_y = self.classifier(z)
        return w_y

    @auto_move_data
    def classification_loss(self, labelled_dataset):  # noqa: D102
        x = labelled_dataset[REGISTRY_KEYS.X_KEY]
        y = labelled_dataset[REGISTRY_KEYS.LABELS_KEY]
        batch_idx = labelled_dataset[REGISTRY_KEYS.BATCH_KEY]
        cont_key = REGISTRY_KEYS.CONT_COVS_KEY
        cont_covs = (
            labelled_dataset[cont_key] if cont_key in labelled_dataset.keys() else None
        )

        cat_key = REGISTRY_KEYS.CAT_COVS_KEY
        cat_covs = (
            labelled_dataset[cat_key] if cat_key in labelled_dataset.keys() else None
        )
        classification_loss = F.cross_entropy(
            self.classify(
                x, batch_index=batch_idx, cat_covs=cat_covs, cont_covs=cont_covs
            ),
            y.view(-1).long(),
        )
        return classification_loss

    def loss(
        self,
        tensors,
        inference_outputs,
        generative_ouputs,
        feed_labels=False,
        kl_weight=1,
        labelled_tensors=None,
        classification_ratio=None,
    ):
        """Compute the loss."""
        px = generative_ouputs["px"]
        qz1 = inference_outputs["qz"]
        z1 = inference_outputs["z"]
        x = tensors[REGISTRY_KEYS.X_KEY]
        batch_index = tensors[REGISTRY_KEYS.BATCH_KEY]

        if feed_labels:
            y = tensors[REGISTRY_KEYS.LABELS_KEY]
        else:
            y = None
        is_labelled = False if y is None else True

        # Enumerate choices of label
        ys, z1s = broadcast_labels(y, z1, n_broadcast=self.n_labels)
        qz2, z2 = self.encoder_z2_z1(z1s, ys)
        pz1_m, pz1_v = self.decoder_z1_z2(z2, ys)
        reconst_loss = -px.log_prob(x).sum(-1)

        # KL Divergence
        mean = torch.zeros_like(qz2.loc)
        scale = torch.ones_like(qz2.scale)

        kl_divergence_z2 = kl(qz2, Normal(mean, scale)).sum(dim=1)
        loss_z1_unweight = -Normal(pz1_m, torch.sqrt(pz1_v)).log_prob(z1s).sum(dim=-1)
        loss_z1_weight = qz1.log_prob(z1).sum(dim=-1)
        if not self.use_observed_lib_size:
            ql = inference_outputs["ql"]
            (
                local_library_log_means,
                local_library_log_vars,
            ) = self._compute_local_library_params(batch_index)

            kl_divergence_l = kl(
                ql,
                Normal(local_library_log_means, torch.sqrt(local_library_log_vars)),
            ).sum(dim=1)
        else:
            kl_divergence_l = 0.0

        if is_labelled:
            loss = reconst_loss + loss_z1_weight + loss_z1_unweight
            kl_locals = {
                "kl_divergence_z2": kl_divergence_z2,
                "kl_divergence_l": kl_divergence_l,
            }
            if labelled_tensors is not None:
                classifier_loss = self.classification_loss(labelled_tensors)
                loss += classifier_loss * classification_ratio
                return LossOutput(
                    loss=loss,
                    reconstruction_loss=reconst_loss,
                    kl_local=kl_locals,
                    extra_metrics={
                        "classification_loss": classifier_loss,
                        "n_labelled_tensors": labelled_tensors[
                            REGISTRY_KEYS.X_KEY
                        ].shape[0],
                    },
                )
            return LossOutput(
                loss=loss,
                reconstruction_loss=reconst_loss,
                kl_local=kl_locals,
            )

        probs = self.classifier(z1)
        reconst_loss += loss_z1_weight + (
            (loss_z1_unweight).view(self.n_labels, -1).t() * probs
        ).sum(dim=1)

        kl_divergence = (kl_divergence_z2.view(self.n_labels, -1).t() * probs).sum(
            dim=1
        )
        kl_divergence += kl(
            Categorical(probs=probs),
            Categorical(probs=self.y_prior.repeat(probs.size(0), 1)),
        )
        kl_divergence += kl_divergence_l

        loss = torch.mean(reconst_loss + kl_divergence * kl_weight)

        if labelled_tensors is not None:
            classifier_loss = self.classification_loss(labelled_tensors)
            loss += classifier_loss * classification_ratio
            return LossOutput(
                loss=loss,
                reconstruction_loss=reconst_loss,
                kl_local=kl_divergence,
                extra_metrics={"classification_loss": classifier_loss},
            )
        return LossOutput(
            loss=loss, reconstruction_loss=reconst_loss, kl_local=kl_divergence
        )
