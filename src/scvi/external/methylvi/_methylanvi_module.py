"""PyTorch module for methylVI for single cell methylation data."""

from collections.abc import Iterable, Sequence
from typing import Literal

import numpy as np
import torch
from torch.distributions import Categorical, Normal
from torch.distributions import kl_divergence as kl
from torch.nn import functional as F

from scvi import REGISTRY_KEYS
from scvi.external.methylvi._base_components import BSSeqModuleMixin
from scvi.external.methylvi._methylvi_module import METHYLVAE
from scvi.module._classifier import Classifier
from scvi.module._utils import broadcast_labels
from scvi.module.base import LossOutput, SupervisedModuleClass, auto_move_data
from scvi.nn import Decoder, Encoder


class METHYLANVAE(SupervisedModuleClass, METHYLVAE, BSSeqModuleMixin):
    """Methylation annotation using variational inference.

    This is an implementation of the MethylANVI model described in :cite:p:`Weinberger2023a`.

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
    likelihood
        One of
        * ``'betabinomial'`` - BetaBinomial distribution
        * ``'binomial'`` - Binomial distribution
    dispersion
        One of the following
        * ``'region'`` - dispersion parameter of BetaBinomial is constant per region across cells
        * ``'region-cell'`` - dispersion can differ for every region in every cell
    y_prior
        If None, initialized to uniform probability over cell types
    labels_groups
        Label group designations
    use_labels_groups
        Whether to use the label groups
    linear_classifier
        If `True`, uses a single linear layer for classification instead of a
        multi-layer perceptron.
    classifier_parameters
        Keyword arguments passed into :class:`~scvi.module.Classifier`.
    use_batch_norm
        Whether to use batch norm in layers
    use_layer_norm
        Whether to use layer norm in layers
    linear_classifier
        If ``True``, uses a single linear layer for classification instead of a
        multi-layer perceptron.
    **model_kwargs
        Keyword args for :class:`~scvi.external.methylvi.METHYLVAE`
    """

    def __init__(
        self,
        n_input: int,
        contexts: Iterable[str],
        num_features_per_context: Iterable[int],
        n_batch: int = 0,
        n_cats_per_cov: Iterable[int] | None = None,
        n_labels: int = 0,
        n_hidden: int = 128,
        n_latent: int = 10,
        n_layers: int = 1,
        dropout_rate: float = 0.1,
        likelihood: Literal["betabinomial", "binomial"] = "betabinomial",
        dispersion: Literal["region", "region-cell"] = "region",
        y_prior=None,
        labels_groups: Sequence[int] = None,
        use_labels_groups: bool = False,
        linear_classifier: bool = False,
        classifier_parameters: dict | None = None,
        use_batch_norm: Literal["encoder", "decoder", "none", "both"] = "both",
        use_layer_norm: Literal["encoder", "decoder", "none", "both"] = "none",
        **model_kwargs,
    ):
        super().__init__(
            n_input=n_input,
            n_hidden=n_hidden,
            n_latent=n_latent,
            n_layers=n_layers,
            n_batch=n_batch,
            n_cats_per_cov=n_cats_per_cov,
            contexts=contexts,
            num_features_per_context=num_features_per_context,
            likelihood=likelihood,
            dispersion=dispersion,
            **model_kwargs,
        )

        classifier_parameters = classifier_parameters or {}
        use_batch_norm_encoder = use_batch_norm == "encoder" or use_batch_norm == "both"
        use_batch_norm_decoder = use_batch_norm == "decoder" or use_batch_norm == "both"
        use_layer_norm_encoder = use_layer_norm == "encoder" or use_layer_norm == "both"
        use_layer_norm_decoder = use_layer_norm == "decoder" or use_layer_norm == "both"

        self.n_labels = n_labels
        # Classifier takes n_latent as input
        cls_parameters = {
            "n_layers": 0 if linear_classifier else n_layers,
            "n_hidden": 0 if linear_classifier else n_hidden,
            "dropout_rate": dropout_rate,
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
    def classify(
        self,
        mc: torch.Tensor,
        cov: torch.Tensor,
        batch_index: torch.Tensor | None = None,
        cont_covs=None,
        cat_covs=None,
        use_posterior_mean: bool = True,
    ) -> torch.Tensor:
        """Forward pass through the encoder and classifier of methylANVI."""
        # get variational parameters via the encoder networks
        # we input both the methylated reads (mc) and coverage (cov)
        return super().classify(
            x=torch.cat((mc, cov), dim=-1),
            batch_index=batch_index,
            cont_covs=cont_covs,
            cat_covs=cat_covs,
            use_posterior_mean=use_posterior_mean,
        )

    def loss(
        self,
        tensors,
        inference_outputs,
        generative_outputs,
        feed_labels=False,
        kl_weight=1,
        labelled_tensors=None,
        classification_ratio=None,
    ):
        """Compute the loss."""
        qz1 = inference_outputs["qz"]
        z1 = inference_outputs["z"]

        if feed_labels:
            y = tensors[REGISTRY_KEYS.LABELS_KEY]
        else:
            y = None
        is_labelled = False if y is None else True

        # Enumerate choices of label
        ys, z1s = broadcast_labels(z1, n_broadcast=self.n_labels)
        qz2, z2 = self.encoder_z2_z1(z1s, ys)
        pz1_m, pz1_v = self.decoder_z1_z2(z2, ys)

        minibatch_size = qz1.loc.size()[0]
        reconst_loss = self._compute_minibatch_reconstruction_loss(
            minibatch_size=minibatch_size,
            tensors=tensors,
            generative_outputs=generative_outputs,
        )

        # KL Divergence
        mean = torch.zeros_like(qz2.loc)
        scale = torch.ones_like(qz2.scale)

        kl_divergence_z2 = kl(qz2, Normal(mean, scale)).sum(dim=1)
        loss_z1_unweight = -Normal(pz1_m, torch.sqrt(pz1_v)).log_prob(z1s).sum(dim=-1)
        loss_z1_weight = qz1.log_prob(z1).sum(dim=-1)

        if is_labelled:
            loss = reconst_loss + loss_z1_weight + loss_z1_unweight
            kl_locals = {
                "kl_divergence_z2": kl_divergence_z2,
            }
            if labelled_tensors is not None:
                ce_loss, true_labels, logits = self.classification_loss(labelled_tensors)
                loss += ce_loss * classification_ratio
                return LossOutput(
                    loss=loss,
                    reconstruction_loss=reconst_loss,
                    kl_local=kl_locals,
                    classification_loss=ce_loss,
                    true_labels=true_labels,
                    logits=logits,
                    extra_metrics={
                        "n_labelled_tensors": labelled_tensors[REGISTRY_KEYS.X_KEY].shape[0],
                    },
                )
            return LossOutput(
                loss=loss,
                reconstruction_loss=reconst_loss,
                kl_local=kl_locals,
            )

        probs = F.softmax(self.classifier(z1), dim=-1)

        reconst_loss += loss_z1_weight + (
            (loss_z1_unweight).view(self.n_labels, -1).t() * probs
        ).sum(dim=1)

        kl_divergence = (kl_divergence_z2.view(self.n_labels, -1).t() * probs).sum(dim=1)
        kl_divergence += kl(
            Categorical(probs=probs),
            Categorical(probs=self.y_prior.repeat(probs.size(0), 1)),
        )

        loss = torch.mean(reconst_loss + kl_divergence * kl_weight)

        if labelled_tensors is not None:
            ce_loss, true_labels, logits = self.classification_loss(labelled_tensors)

            loss += ce_loss * classification_ratio
            return LossOutput(
                loss=loss,
                reconstruction_loss=reconst_loss,
                kl_local=kl_divergence,
                classification_loss=ce_loss,
                true_labels=true_labels,
                logits=logits,
            )
        return LossOutput(loss=loss, reconstruction_loss=reconst_loss, kl_local=kl_divergence)
