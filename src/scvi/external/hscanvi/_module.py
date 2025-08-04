from __future__ import annotations

from typing import TYPE_CHECKING

import numpy as np
import torch
from torch.distributions import Categorical, Normal
from torch.distributions import kl_divergence as kl
from torch.nn import functional as F

from scvi import REGISTRY_KEYS
from scvi.module.base import LossOutput, SupervisedModuleClass
from scvi.nn import Decoder, Encoder

from scvi.module.base import auto_move_data

from scvi.module._classifier import Classifier
from scvi.module._utils import broadcast_labels
from scvi.module._scanvae import SCANVAE

if TYPE_CHECKING:
    from collections.abc import Iterable, Sequence
    from typing import Literal

    from torch.distributions import Distribution


class HSCANVAE(SCANVAE):

    def __init__(
        self,
        n_labels_multilabel: int = 0,
        **scanvae_kwargs,
    ):
        super().__init__(**scanvae_kwargs)

        classifier_parameters = {}

        self.n_labels_multilabel = n_labels_multilabel

        linear_classifier = scanvae_kwargs.get("linear_classifier", False)
        n_layers = scanvae_kwargs.get("n_layers", 1)
        n_hidden = scanvae_kwargs.get("n_hidden", 128)
        dropout_rate = scanvae_kwargs.get("dropout_rate", 0.1)
        use_batch_norm_encoder = (
            scanvae_kwargs.get("use_batch_norm", "both") == "encoder"
            or scanvae_kwargs.get("use_batch_norm", "both") == "both"
        )
        use_layer_norm_encoder = (
            scanvae_kwargs.get("use_layer_norm", "both") == "encoder"
            or scanvae_kwargs.get("use_layer_norm", "both") == "both"
        )
        n_latent = scanvae_kwargs.get("n_latent", 10)

        cls_parameters = {
            "n_layers": 0 if [linear_classifier] else n_layers,
            "n_hidden": 0 if linear_classifier else n_hidden,
            "dropout_rate": dropout_rate,
            "logits": True,
        }
        cls_parameters.update(classifier_parameters)
        self.classifier = Classifier(
            n_latent,
            n_labels=n_labels_multilabel,
            use_batch_norm=use_batch_norm_encoder,
            use_layer_norm=use_layer_norm_encoder,
            **cls_parameters,
        )

    @auto_move_data
    def classification_loss(
        self, labelled_dataset: dict[str, torch.Tensor]
    ) -> tuple[torch.Tensor, torch.Tensor, torch.Tensor]:
        inference_inputs = self._get_inference_input(
            labelled_dataset
        )  # (n_obs, n_vars)
        data_inputs = {
            key: inference_inputs[key]
            for key in inference_inputs.keys()
            if key not in ["batch_index", "cont_covs", "cat_covs", "panel_index"]
        }

        y = labelled_dataset[REGISTRY_KEYS.LABELS_KEY]  # (n_obs, 1)
        batch_idx = labelled_dataset[REGISTRY_KEYS.BATCH_KEY]
        cont_key = REGISTRY_KEYS.CONT_COVS_KEY
        cont_covs = (
            labelled_dataset[cont_key] if cont_key in labelled_dataset.keys() else None
        )

        cat_key = REGISTRY_KEYS.CAT_COVS_KEY
        cat_covs = (
            labelled_dataset[cat_key] if cat_key in labelled_dataset.keys() else None
        )
        # NOTE: prior to v1.1, this method returned probabilities per label by
        # default, see #2301 for more details
        logits = self.classify(
            **data_inputs, batch_index=batch_idx, cat_covs=cat_covs, cont_covs=cont_covs
        )  # (n_obs, n_labels)
        ce_loss = F.cross_entropy(
            logits,
            y.view(-1).long(),
        )
        return ce_loss, y, logits

    def loss(
        self,
        tensors: dict[str, torch.Tensor],
        inference_outputs: dict[str, torch.Tensor | Distribution | None],
        generative_ouputs: dict[str, Distribution | None],
        kl_weight: float = 1.0,
        labelled_tensors: dict[str, torch.Tensor] | None = None,
        classification_ratio: float | None = None,
    ):
        """Compute the loss."""
        px: Distribution = generative_ouputs["px"]
        qz1: torch.Tensor = inference_outputs["qz"]
        z1: torch.Tensor = inference_outputs["z"]
        x: torch.Tensor = tensors[REGISTRY_KEYS.X_KEY]
        batch_index: torch.Tensor = tensors[REGISTRY_KEYS.BATCH_KEY]

        ys, z1s = broadcast_labels(z1, n_broadcast=self.n_labels)
        qz2, z2 = self.encoder_z2_z1(z1s, ys)
        pz1_m, pz1_v = self.decoder_z1_z2(z2, ys)
        reconst_loss = -px.log_prob(x).sum(-1)

        # KL Divergence
        mean = torch.zeros_like(qz2.loc)
        scale = torch.ones_like(qz2.scale)

        kl_divergence_z2 = kl(qz2, Normal(mean, scale)).sum(dim=-1)
        loss_z1_unweight = -Normal(pz1_m, torch.sqrt(pz1_v)).log_prob(z1s).sum(dim=-1)
        loss_z1_weight = qz1.log_prob(z1).sum(dim=-1)

        probs = self.classifier(z1)
        if self.classifier.logits:
            probs = F.softmax(probs, dim=-1)

        if z1.ndim == 2:
            loss_z1_unweight_ = loss_z1_unweight.view(self.n_labels, -1).t()
            kl_divergence_z2_ = kl_divergence_z2.view(self.n_labels, -1).t()
        else:
            loss_z1_unweight_ = torch.transpose(
                loss_z1_unweight.view(z1.shape[0], self.n_labels, -1), -1, -2
            )
            kl_divergence_z2_ = torch.transpose(
                kl_divergence_z2.view(z1.shape[0], self.n_labels, -1), -1, -2
            )
        reconst_loss += loss_z1_weight + (loss_z1_unweight_ * probs).sum(dim=-1)
        kl_divergence = (kl_divergence_z2_ * probs).sum(dim=-1)
        kl_divergence += kl(
            Categorical(probs=probs),
            Categorical(
                probs=(
                    self.y_prior.repeat(probs.size(0), probs.size(1), 1)
                    if len(probs.size()) == 3
                    else self.y_prior.repeat(probs.size(0), 1)
                )
            ),
        )

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
            kl_divergence_l = torch.zeros_like(kl_divergence)

        kl_divergence += kl_divergence_l

        loss = torch.mean(reconst_loss + kl_divergence * kl_weight)

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
                reconstruction_loss=reconst_loss,
                kl_local=kl_divergence,
                classification_loss=ce_loss,
                true_labels=true_labels,
                logits=logits,
                extra_metrics=extra_metrics_payload,
            )
        return LossOutput(
            loss=loss,
            reconstruction_loss=reconst_loss,
            kl_local=kl_divergence,
            extra_metrics=extra_metrics_payload,
        )
