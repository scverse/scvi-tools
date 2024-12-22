from typing import Literal

import torch
from torch.distributions import Categorical, Independent, MixtureSameFamily, Normal
from torch.distributions import kl_divergence as kl
from torch.nn import functional as F

from scvi import REGISTRY_KEYS
from scvi.module import SCANVAE
from scvi.module._utils import broadcast_labels
from scvi.module.base import LossOutput, auto_move_data

from ._base_components import HierarchicalLossNetwork, MultiBatchClassifier


class MUANVAE(SCANVAE):
    """
    Single-cell multiple-annotation using variational inference.

    This is a re-implementation of a hierarchical cell-type annotation model
    inspired from scANVI model described in [Xu21]_,.

    Parameters
    ----------
    n_input
        Number of input genes
    n_batch
        Number of batches
    n_fine_labels
        Number of fine labels
    num_classes
        Number of labels per class organized in a hierarchical list
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
    conditioning_class
        index of the class conditioning the second latent space (0 being the coarse class, 1 being the fine class). Default : coarse labels.
    hierarchy_matrix
        Matrix representing the hierarchy, computed by scATVI. If None, no hierarchical y prior update in the loss.
    use_batch_norm
        Whether to use batch norm in layers
    use_layer_norm
        Whether to use layer norm in layers
    prior_z1
        Whether to use MoG or simple Gaussian for prior of z1
    **vae_kwargs
        Keyword args for :class:`~scvi.module.VAE`
    """

    def __init__(
        self,
        n_input: int,
        num_classes: list,
        hierarchy_dict: dict,
        n_batch: int = 0,
        n_site: int = 0,
        n_fine_labels: int = 0,
        n_hidden: int = 128,
        n_latent: int = 10,
        n_layers: int = 1,
        n_continuous_cov: int = 0,
        n_cats_per_cov: list[int] | None = None,
        dropout_rate: float = 0.1,
        dispersion: str = "gene",
        log_variational: bool = True,
        gene_likelihood: str = "nb",
        classifier_parameters: dict = dict(),
        classifier_parameters_muanvae: dict = dict(),
        use_batch_norm: Literal["encoder", "decoder", "none", "both"] = "none",
        use_layer_norm: Literal["encoder", "decoder", "none", "both"] = "both",
        conditioning_class: int = -1,
        mog_class: int = 0,
        hierarchy_matrix=None,
        update_yprior=True,
        prior_z1: str = "gaussian",
        eps_yprior: float = 1e-6,
        **scanvae_kwargs,
    ):
        self.conditioning_class = conditioning_class
        self.mog_class = mog_class
        self.site_specific_classifier = n_site > 1
        self.n_site = n_site
        self.num_classes = num_classes
        self.update_yprior = update_yprior
        self.n_labels_conditioning = num_classes[self.conditioning_class]
        self.n_fine_labels = n_fine_labels
        self.hiearchy_dict = hierarchy_dict

        classifier_parameters = classifier_parameters or {}

        cls_parameters = {
            "n_layers": n_layers,
            "n_hidden": n_hidden,
            "dropout_rate": dropout_rate,
            "logits": True,
        }
        cls_parameters.update(classifier_parameters)

        super().__init__(
            n_input,
            n_batch=n_batch,
            n_labels=self.n_fine_labels,
            n_hidden=n_hidden,
            n_latent=n_latent,
            n_layers=n_layers,
            n_continuous_cov=n_continuous_cov,
            n_cats_per_cov=n_cats_per_cov,
            dropout_rate=dropout_rate,
            dispersion=dispersion,
            log_variational=log_variational,
            gene_likelihood=gene_likelihood,
            classifier_parameters=classifier_parameters,
            use_batch_norm=use_batch_norm,
            use_layer_norm=use_layer_norm,
            **scanvae_kwargs,
        )
        cls_parameters.update(classifier_parameters_muanvae)
        self.num_classes = num_classes
        self.total_level = len(self.num_classes)  # depth of hierarchy
        self.prior_z1 = prior_z1

        self.classifier = HierarchicalLossNetwork(
            n_input=self.n_latent,
            num_classes=self.num_classes[:-1],
            **cls_parameters,
        )
        # the site-specific classifiers must have layer norm (batch norm does not work if there is 1 single observation from a batch in a minibatch)
        self.multi_classifier_fine = MultiBatchClassifier(
            n_input=self.n_latent,
            n_sites=n_site,
            n_labels=n_fine_labels,
            use_batch_norm=False,
            use_layer_norm=True,
            **cls_parameters,
        )

        # register y_prior on the fine labels
        if not self.update_yprior:
            hierarchy_matrix = [
                torch.tensor(hierarchy_matrix[i].sum(0) > 0, dtype=torch.float)
                for i in range(len(num_classes) - 1)
            ]
        for i in range(0, self.total_level):
            if i == self.total_level - 1:
                self.y_prior_fine = torch.nn.ParameterList(
                    [
                        torch.nn.Parameter(  #
                            hierarchy_matrix[i - 1][:, :, site] + eps_yprior,
                            requires_grad=False,
                        )
                        for site in range(n_site)
                    ]
                )
            elif i == 0:
                self.register_buffer(
                    f"y_prior_{i}",
                    torch.nn.Parameter(  #
                        torch.full([num_classes[0]], 1 / num_classes[0], dtype=torch.float),
                        requires_grad=False,
                    ),
                )
            else:
                self.register_buffer(
                    f"y_prior_{i}",
                    torch.nn.Parameter(  #
                        torch.tensor(hierarchy_matrix[i - 1] + eps_yprior, dtype=torch.float),
                        requires_grad=False,
                    ),
                )
        if self.prior_z1 == "mog":
            self.register_parameter(
                "prior_z1_means",
                torch.nn.Parameter(torch.randn([num_classes[self.mog_class], n_latent])),
            )
            self.register_parameter(
                "prior_z1_scales",
                torch.nn.Parameter(torch.zeros([num_classes[self.mog_class], n_latent])),
            )
            self.register_parameter(
                "prior_z1_logits", torch.nn.Parameter(torch.ones([num_classes[self.mog_class]]))
            )
        if self.prior_z1 == "mog_celltype":
            self.register_parameter(
                "prior_z1_means",
                torch.nn.Parameter(torch.zeros([1, num_classes[self.mog_class], n_latent])),
            )
            self.register_parameter(
                "prior_z1_scales",
                torch.nn.Parameter(torch.zeros([1, num_classes[self.mog_class], n_latent])),
            )

    @auto_move_data
    def classify(
        self,
        x,
        batch_index=None,
        site_index=None,
        cont_covs=None,
        cat_covs=None,
        site_to_predict: int | None = None,
        precomputed_z: torch.Tensor | None = None,
    ):
        """
        Classify cells using the model.

        Parameters
        ----------
        site_to_predict
            For cross prediction purposes : index of the fine classifier to use when cross-predicting labels.
            If None, normal prediction occurs, which is the case during training.
            If not None, cross classification on only fine layer with this batch-specific classifier.
        precomputed_z
            Precomputed z1 latent space. If None, z1 is computed from x.
        """
        if precomputed_z is not None:
            z = precomputed_z
        else:
            if self.log_variational:  # for numerical stability
                x = torch.log(1 + x)

            if cont_covs is not None and self.encode_covariates:
                encoder_input = torch.cat((x, cont_covs), dim=-1)
            else:
                encoder_input = x
            if cat_covs is not None and self.encode_covariates:
                categorical_input = torch.split(cat_covs, 1, dim=1)
            else:
                categorical_input = ()
            qz, _ = self.z_encoder(
                encoder_input, batch_index, *categorical_input
            )  # q(z1|x)  without the var   qz_v
            z = qz.rsample()
        if site_to_predict is not None:
            site_to_predict_index = torch.full((z.shape[0], 1), site_to_predict, dtype=torch.int32)
            probs_fine, _ = self.multi_classifier_fine(z, site_to_predict_index)
            return probs_fine, _
        probs_classifier, logits_classifier = self.classifier(z)
        probs_fine, logits_fine = self.multi_classifier_fine(z, site_index)
        probs_classifier += [probs_fine]
        logits_classifier += [logits_fine]

        return probs_classifier, logits_classifier

    @auto_move_data
    def classification(
        self,
        tensors,
        return_classifier_loss=False,
        site_to_predict=None,
        precomputed_z=None,
    ):
        x = tensors[REGISTRY_KEYS.X_KEY]
        y = tensors["label_hierarchy"]
        batch_idx = tensors[REGISTRY_KEYS.BATCH_KEY]
        site_idx = tensors[REGISTRY_KEYS.SITE_KEY]
        cont_covs = tensors.get(REGISTRY_KEYS.CONT_COVS_KEY, None)
        cat_covs = tensors.get(REGISTRY_KEYS.CAT_COVS_KEY, None)

        probs, logits = self.classify(
            x,
            batch_index=batch_idx,
            site_index=site_idx,
            cat_covs=cat_covs,
            cont_covs=cont_covs,
            site_to_predict=site_to_predict,
            precomputed_z=precomputed_z,
        )
        if not return_classifier_loss:
            return probs, logits

        classification_loss = 0
        for l in range(self.total_level):
            labels_curr_level = y[:, l].view(-1).long()
            classification_loss += F.cross_entropy(
                logits[l],
                labels_curr_level,
                ignore_index=self.num_classes[l],
                reduction="mean"
            )
        true_labels_fine = torch.unsqueeze(labels_curr_level, 1)
        return classification_loss, true_labels_fine, logits[-1]

    def loss(
        self,
        tensors,
        inference_outputs,
        generative_ouputs,
        kl_weight=1,
        labelled_tensors=None,
        classification_ratio=None,
        bg_classifier_ratio=0.0,
        loss_z1_factor=0.1,
        weighting_mog=1.0,
        warmup_model=False,  # if true, trains a model without cell-type classification first.
    ):
        """Compute the loss."""
        px = generative_ouputs["px"]
        qz1 = inference_outputs["qz"]
        z1 = inference_outputs["z"]
        x = tensors[REGISTRY_KEYS.X_KEY]
        y = tensors[REGISTRY_KEYS.LABELS_KEY]
        site_index = tensors[REGISTRY_KEYS.SITE_KEY]

        is_labelled = False if y is None else True

        # Enumerate choices of label
        ys, z1s = broadcast_labels(z1, n_broadcast=self.n_labels_conditioning)
        qz2, z2 = self.encoder_z2_z1(z1s, ys)
        pz1_m, pz1_v = self.decoder_z1_z2(z2, ys)
        reconst_loss = -px.log_prob(x).sum(-1)

        # KL Divergence
        mean = torch.zeros_like(qz2.loc)
        scale = torch.ones_like(qz2.scale)

        kl_divergence_z2 = kl(qz2, Normal(mean, scale)).sum(dim=1)
        loss_z1_unweight = -Normal(pz1_m, torch.sqrt(pz1_v)).log_prob(z1s).sum(dim=-1)
        loss_z1_weight = qz1.log_prob(z1).sum(dim=-1)

        probs, logits = self.classification(tensors, precomputed_z=z1)
        probs_conditioning, logits_conditioning = (
            probs[self.conditioning_class],
            logits[self.conditioning_class],
        )

        if not warmup_model:
            reconst_loss += (
                (
                    loss_z1_weight
                    + ((loss_z1_unweight).view(self.n_labels_conditioning, -1).t()
                    * probs_conditioning).sum(dim=1)
                )
                * kl_weight
                * loss_z1_factor
            )

        if self.prior_z1 == "mog":
            cats = Categorical(logits=self.prior_logits)
            normal_dists = Independent(
                Normal(self.prior_means, torch.exp(self.prior_log_scales) + 1e-4),
                1,
            )
            prior = MixtureSameFamily(cats, normal_dists)
            u = qz1.rsample(sample_shape=(30,))
            # (sample, n_obs, n_latent) -> (sample, n_obs,)
            kl_divergence = -(prior.log_prob(u) - qz1.log_prob(u).sum(-1)).mean(0)
        elif self.prior_z1 == "mog_celltype":
            if warmup_model:
                # Assigns zero meaning equal weight to all unlabeled cells. Otherwise biases to sample from respective MoG.
                logits_input = torch.nn.functional.one_hot(
                    y[:, self.mog_class].ravel().long(), self.num_classes[self.mog_class] + 1
                ).float()[:, :-1]
                cats = Categorical(logits=10 * logits_input)
            else:
                cats = Categorical(logits=logits_conditioning)
            normal_dists = torch.distributions.Independent(
                Normal(
                    self.prior_z1_means.expand(x.shape[0], -1, -1),
                    torch.exp(self.prior_z1_scales).expand(x.shape[0], -1, -1) + 1e-2,
                ),
                reinterpreted_batch_ndims=1,
            )

            prior = MixtureSameFamily(cats, normal_dists)
            u = qz1.rsample(sample_shape=(30,))
            # (sample, n_obs, n_latent) -> (sample, n_obs,)
            kl_z = -(prior.log_prob(u) - qz1.log_prob(u).sum(-1)).mean(0)
            kl_divergence = weighting_mog * kl_z
        else:
            prior = Normal(torch.zeros_like(qz1.loc), torch.ones_like(qz1.loc))
            kl_z = 0
            kl_divergence = 0

        probs_prior, _ = self.classification(tensors, precomputed_z=prior.sample())
        kl_divergence_cat = 0

        if not warmup_model:
            for i in range(0, self.total_level):
                y_prior_ = (
                    torch.stack(
                        [self.y_prior_fine[idx] for idx in site_index.ravel().long()], dim=0
                    )
                    if i == self.total_level - 1
                    else getattr(self, f"y_prior_{i}")
                )
                if self.update_yprior and i > 0:
                    if i == self.total_level - 1:
                        # Shape batch, coarse; batch, coarse, finer -> batch, finer
                        y_prior_ = torch.einsum("bc,bcf->bf", probs[i - 1], y_prior_)
                    else:
                        y_prior_ = torch.einsum("bc,cf->bf", probs[i - 1], y_prior_)

                kl_divergence_cat += kl(
                    Categorical(probs=probs[i]),
                    Categorical(probs=y_prior_),
                )
            for i in range(0, self.total_level):
                if i == self.total_level - 1:
                    y_prior_ = torch.stack(
                        [self.y_prior_fine[i] for i in site_index.ravel().long()], dim=0
                    )
                else:
                    y_prior_ = self.__getattr__("y_prior_" + str(i))
                if self.update_yprior and i > 0:
                    if i == self.total_level - 1:
                        # Shape batch, coarse; batch, coarse, finer -> batch, finer
                        y_prior_ = torch.einsum("bc,bcf->bf", probs_prior[i - 1], y_prior_)
                    else:
                        y_prior_ = torch.einsum("bc,cf->bf", probs_prior[i - 1], y_prior_)

                kl_divergence_cat += bg_classifier_ratio * kl(
                    Categorical(probs=probs_prior[i]),
                    Categorical(probs=y_prior_),
                )

        kl_divergence += kl_divergence_cat

        loss = torch.mean(reconst_loss + kl_divergence * kl_weight)

        if labelled_tensors is not None:
            # We filter cells with unlabeled_category in loss.
            ce_loss, fine_true_labels, logits_fine = self.classification(
                tensors, return_classifier_loss=True
            )
            if not warmup_model:
                loss += ce_loss * classification_ratio
            return LossOutput(
                loss=loss,
                reconstruction_loss=reconst_loss,
                kl_local=kl_divergence,
                classification_loss=ce_loss,
                true_labels=fine_true_labels,
                logits=logits_fine,
            )
        return LossOutput(loss=loss, reconstruction_loss=reconst_loss, kl_local=kl_divergence)
