"""PyTorch module for methylVI for single cell methylation data."""

from collections.abc import Iterable, Sequence
from typing import Literal

import numpy as np
import torch
import torch.nn as nn
from torch.distributions import Binomial, Categorical, Normal
from torch.distributions import kl_divergence as kl
from torch.nn import functional as F

from scvi import REGISTRY_KEYS
from scvi.distributions import BetaBinomial
from scvi.external.methylvi import METHYLVI_REGISTRY_KEYS, DecoderMETHYLVI
from scvi.external.methylvi._utils import _context_cov_key, _context_mc_key
from scvi.module._classifier import Classifier
from scvi.module._utils import broadcast_labels
from scvi.module.base import BaseModuleClass, LossOutput, auto_move_data
from scvi.nn import Decoder, Encoder

TensorDict = dict[str, torch.Tensor]


class METHYLVAE(BaseModuleClass):
    """PyTorch module for methylVI.

    Parameters
    ----------
    n_input
        Total number of input genomic regions
    contexts
        List of methylation contexts (e.g. ["mCG", "mCH"])
    num_features_per_context
        Number of features corresponding to each context
    n_batch
        Number of batches, if 0, no batch correction is performed
    n_cats_per_cov
        Number of categories for each extra categorical covariate
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
        * ``'region-cell'`` - dispersion can differ for every region in every cell
    """

    def __init__(
        self,
        n_input: int,
        contexts: Iterable[str],
        num_features_per_context: Iterable[int],
        n_batch: int = 0,
        n_cats_per_cov: Iterable[int] | None = None,
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
        self.contexts = contexts
        self.num_features_per_context = num_features_per_context

        cat_list = [n_batch] + list([] if n_cats_per_cov is None else n_cats_per_cov)

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
        for context, num_features in zip(contexts, num_features_per_context, strict=False):
            self.decoders[context] = DecoderMETHYLVI(
                n_latent,
                num_features,
                n_cat_list=cat_list,
                n_layers=n_layers,
                n_hidden=n_hidden,
            )

        if self.dispersion == "region":
            self.px_gamma = torch.nn.ParameterDict(
                {
                    context: nn.Parameter(torch.randn(num_features))
                    for (context, num_features) in zip(
                        contexts, num_features_per_context, strict=False
                    )
                }
            )

    def _get_inference_input(self, tensors):
        """Parse the dictionary to get appropriate args"""
        mc = torch.cat(
            [tensors[_context_mc_key(context)] for context in self.contexts],
            dim=1,
        )
        cov = torch.cat(
            [tensors[_context_cov_key(context)] for context in self.contexts],
            dim=1,
        )

        batch_index = tensors[REGISTRY_KEYS.BATCH_KEY]
        cat_key = REGISTRY_KEYS.CAT_COVS_KEY
        cat_covs = tensors[cat_key] if cat_key in tensors.keys() else None

        input_dict = {
            METHYLVI_REGISTRY_KEYS.MC_KEY: mc,
            METHYLVI_REGISTRY_KEYS.COV_KEY: cov,
            "batch_index": batch_index,
            "cat_covs": cat_covs,
        }
        return input_dict

    def _get_generative_input(self, tensors, inference_outputs):
        z = inference_outputs["z"]
        batch_index = tensors[REGISTRY_KEYS.BATCH_KEY]
        cat_key = REGISTRY_KEYS.CAT_COVS_KEY
        cat_covs = tensors[cat_key] if cat_key in tensors.keys() else None

        input_dict = {
            "z": z,
            "batch_index": batch_index,
            "cat_covs": cat_covs,
        }
        return input_dict

    @auto_move_data
    def inference(self, mc, cov, batch_index, cat_covs=None, n_samples=1):
        """
        High level inference method.

        Runs the inference (encoder) model.
        """
        # log the inputs to the variational distribution for numerical stability
        mc_ = torch.log(1 + mc)
        cov_ = torch.log(1 + cov)

        # get variational parameters via the encoder networks
        # we input both the methylated reads (mc) and coverage (cov)
        methylation_input = torch.cat((mc_, cov_), dim=-1)
        if cat_covs is not None:
            categorical_input = torch.split(cat_covs, 1, dim=1)
        else:
            categorical_input = ()

        qz, z = self.z_encoder(methylation_input, batch_index, *categorical_input)
        if n_samples > 1:
            z = qz.sample((n_samples,))

        outputs = {"z": z, "qz": qz}
        return outputs

    @auto_move_data
    def generative(self, z, batch_index, cat_covs=None):
        """Runs the generative model."""
        # form the parameters of the BetaBinomial likelihood
        px_mu, px_gamma = {}, {}
        if cat_covs is not None:
            categorical_input = torch.split(cat_covs, 1, dim=1)
        else:
            categorical_input = ()

        for context in self.contexts:
            px_mu[context], px_gamma[context] = self.decoders[context](
                self.dispersion, z, batch_index, *categorical_input
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

        for context in self.contexts:
            px_mu = generative_outputs["px_mu"][context]
            px_gamma = generative_outputs["px_gamma"][context]
            mc = tensors[f"{context}_{METHYLVI_REGISTRY_KEYS.MC_KEY}"]
            cov = tensors[f"{context}_{METHYLVI_REGISTRY_KEYS.COV_KEY}"]

            if self.dispersion == "region":
                px_gamma = torch.sigmoid(self.px_gamma[context])

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
        for context in self.contexts:
            px_mu = generative_outputs["px_mu"][context]
            px_gamma = generative_outputs["px_gamma"][context]
            cov = tensors[f"{context}_{METHYLVI_REGISTRY_KEYS.COV_KEY}"]

            if self.dispersion == "region":
                px_gamma = torch.sigmoid(self.px_gamma[context])

            if self.likelihood == "binomial":
                dist = Binomial(probs=px_mu, total_count=cov)
            elif self.likelihood == "betabinomial":
                dist = BetaBinomial(mu=px_mu, gamma=px_gamma, total_count=cov)

            if n_samples > 1:
                exprs_ = dist.sample()
                exprs[context] = exprs_.permute(
                    [1, 2, 0]
                ).cpu()  # Shape : (n_cells_batch, n_regions, n_samples)
            else:
                exprs[context] = dist.sample().cpu()

        return exprs


class METHYLANVAE(METHYLVAE):
    """Single-cell annotation using variational inference.

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
    likelihood
        One of
        * ``'betabinomial'`` - BetaBinomial distribution
        * ``'binomial'`` - Binomial distribution
    dispersion
        One of the following
        * ``'region'`` - dispersion parameter of BetaBinomial is constant per region across cells
        * ``'region-cell'`` - dispersion can differ for every region in every cell
    log_variational
        Log(data+1) prior to encoding for numerical stability. Not normalization.
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
        cat_covs=None,
        use_posterior_mean: bool = True,
    ) -> torch.Tensor:
        """Forward pass through the encoder and classifier.

        Parameters
        ----------
        x
            Tensor of shape ``(n_obs, n_vars)``.
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
        # log the inputs to the variational distribution for numerical stability
        mc_ = torch.log(1 + mc)
        cov_ = torch.log(1 + cov)

        # get variational parameters via the encoder networks
        # we input both the methylated reads (mc) and coverage (cov)
        methylation_input = torch.cat((mc_, cov_), dim=-1)
        if cat_covs is not None:
            categorical_input = torch.split(cat_covs, 1, dim=1)
        else:
            categorical_input = ()

        qz, z = self.z_encoder(methylation_input, batch_index, *categorical_input)
        z = qz.loc if use_posterior_mean else z

        if self.use_labels_groups:
            w_g = self.classifier_groups(z)
            unw_y = self.classifier(z)
            w_y = torch.zeros_like(unw_y)
            for i, group_index in enumerate(self.groups_index):
                unw_y_g = unw_y[:, group_index]
                w_y[:, group_index] = unw_y_g / (unw_y_g.sum(dim=-1, keepdim=True) + 1e-8)
                w_y[:, group_index] *= w_g[:, [i]]
        else:
            w_y = self.classifier(z)
        return w_y

    @auto_move_data
    def classification_loss(self, labelled_dataset):
        """Computes scANVI-style classification loss."""
        inference_inputs = self._get_inference_input(labelled_dataset)  # (n_obs, n_vars)
        mc = inference_inputs[METHYLVI_REGISTRY_KEYS.MC_KEY]
        cov = inference_inputs[METHYLVI_REGISTRY_KEYS.COV_KEY]
        y = labelled_dataset[REGISTRY_KEYS.LABELS_KEY]  # (n_obs, 1)
        batch_idx = labelled_dataset[REGISTRY_KEYS.BATCH_KEY]
        cat_covs = inference_inputs["cat_covs"]

        logits = self.classify(
            mc,
            cov,
            batch_index=batch_idx,
            cat_covs=cat_covs,
        )  # (n_obs, n_labels)
        ce_loss = F.cross_entropy(
            logits,
            y.view(-1).long(),
        )
        return ce_loss, y, logits

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
        reconst_loss = torch.zeros(minibatch_size).to(self.device)

        for context in self.contexts:
            px_mu = generative_outputs["px_mu"][context]
            px_gamma = generative_outputs["px_gamma"][context]
            mc = tensors[f"{context}_mc"]
            cov = tensors[f"{context}_cov"]

            if self.dispersion == "region":
                px_gamma = torch.sigmoid(self.px_gamma[context])

            if self.likelihood == "binomial":
                dist = Binomial(probs=px_mu, total_count=cov)
            elif self.likelihood == "betabinomial":
                dist = BetaBinomial(mu=px_mu, gamma=px_gamma, total_count=cov)

            reconst_loss += -dist.log_prob(mc).sum(dim=-1)

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

        probs = self.classifier(z1)
        if self.classifier.logits:
            probs = F.softmax(probs, dim=-1)

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
