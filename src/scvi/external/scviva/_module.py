from __future__ import annotations

import logging
from dataclasses import dataclass
from typing import TYPE_CHECKING

import torch
from torch.nn.functional import one_hot

from scvi import REGISTRY_KEYS
from scvi.module import VAE, Classifier
from scvi.module._constants import MODULE_KEYS
from scvi.module.base import (
    LossOutput,
    auto_move_data,
)

from ._components import DirichletDecoder, Encoder, NicheDecoder
from ._constants import SCVIVA_MODULE_KEYS, SCVIVA_REGISTRY_KEYS

if TYPE_CHECKING:
    from typing import Literal

    import numpy as np
    from torch.distributions import Distribution

    from scvi._types import LossRecord

logger = logging.getLogger(__name__)


class nicheVAE(VAE):
    """Variational auto-encoder with niche decoders :cite:p:`Levy25`.

    Parameters
    ----------
    n_input
        Number of input features.
    n_batch
        Number of batches. If ``0``, no batch correction is performed.
    n_labels
        Number of labels.
    n_hidden
        Number of nodes per hidden layer. Passed into :class:`~scvi.nn.Encoder` and
        :class:`~scvi.nn.DecoderSCVI`.
    n_latent
        Dimensionality of the latent space.
    n_layers
        Number of hidden layers. Passed into :class:`~scvi.nn.Encoder` and
        :class:`~scvi.nn.DecoderSCVI`.
    n_layers_niche
        Number of hidden layers in the niche state decoder.
    n_layers_compo
        Number of hidden layers in the composition decoder.
    n_hidden_niche
        Number of nodes per hidden layer in the niche state decoder.
    n_hidden_compo
        Number of nodes per hidden layer in the composition decoder.
    n_continuous_cov
        Number of continuous covariates.
    n_cats_per_cov
        A list of integers containing the number of categories for each categorical covariate.
    dropout_rate
        Dropout rate. Passed into :class:`~scvi.nn.Encoder` but not :class:`~scvi.nn.DecoderSCVI`.
    dispersion
        Flexibility of the dispersion parameter when ``gene_likelihood`` is either ``"nb"`` or
        ``"zinb"``. One of the following:

        * ``"gene"``: parameter is constant per gene across cells.
        * ``"gene-batch"``: parameter is constant per gene per batch.
        * ``"gene-label"``: parameter is constant per gene per label.
        * ``"gene-cell"``: parameter is constant per gene per cell.
    log_variational
        If ``True``, use :func:`~torch.log1p` on input data before encoding for numerical stability
        (not normalization).
    gene_likelihood
        Distribution to use for reconstruction in the generative process. One of the following:

        * ``"nb"``: :class:`~scvi.distributions.NegativeBinomial`.
        * ``"zinb"``: :class:`~scvi.distributions.ZeroInflatedNegativeBinomial`.
        * ``"poisson"``: :class:`~scvi.distributions.Poisson`.
    latent_distribution
        Distribution to use for the latent space. One of the following:

        * ``"normal"``: isotropic normal.
        * ``"ln"``: logistic normal with normal params N(0, 1).
    niche_likelihood
        Distribution to use for the niche state. One of the following:

        * ``"poisson"``: :class:`~torch.distributions.Poisson`.
        * ``"gaussian"``: :class:`~torch.distributions.Normal`.
        Default is ``"gaussian"`` and Poisson should be used if the niche state is count data.
    cell_rec_weight
        Weight of the cell reconstruction loss.
    latent_kl_weight
        Weight of the latent KL divergence.
    spatial_weight
        Weight of the spatial losses
    prior_mixture
        If ``True``, use a mixture of Gaussians for the latent space. Else, use unimodal Gaussian.
    prior_mixture_k
        Number of components in the Gaussian mixture.
    semisupervised
        If ``True``, use a classifier to predict cell type labels from the latent space.
    linear_classifier
        If ``True``, use a linear classifier. Else, use a neural network.
    inpute_covariates_niche_decoder
        If ``True``, covariates are concatenated to the input of the niche state decoder.
    encode_covariates
        If ``True``, covariates are concatenated to gene expression prior to passing through
        the encoder(s). Else, only gene expression is used.
    deeply_inject_covariates
        If ``True`` and ``n_layers > 1``, covariates are concatenated to the outputs of hidden
        layers in the encoder(s) (if ``encoder_covariates`` is ``True``) and the decoder prior to
        passing through the next layer.
    batch_representation
        ``EXPERIMENTAL`` Method for encoding batch information. One of the following:

        * ``"one-hot"``: represent batches with one-hot encodings.
        * ``"embedding"``: represent batches with continuously-valued embeddings using
          :class:`~scvi.nn.Embedding`.

        Note that batch representations are only passed into the encoder(s) if
        ``encode_covariates`` is ``True``.
    use_batch_norm
        Specifies where to use :class:`~torch.nn.BatchNorm1d` in the model. One of the following:

        * ``"none"``: don't use batch norm in either encoder(s) or decoder.
        * ``"encoder"``: use batch norm only in the encoder(s).
        * ``"decoder"``: use batch norm only in the decoder.
        * ``"both"``: use batch norm in both encoder(s) and decoder.

        Note: if ``use_layer_norm`` is also specified, both will be applied (first
        :class:`~torch.nn.BatchNorm1d`, then :class:`~torch.nn.LayerNorm`).
    use_layer_norm
        Specifies where to use :class:`~torch.nn.LayerNorm` in the model. One of the following:

        * ``"none"``: don't use layer norm in either encoder(s) or decoder.
        * ``"encoder"``: use layer norm only in the encoder(s).
        * ``"decoder"``: use layer norm only in the decoder.
        * ``"both"``: use layer norm in both encoder(s) and decoder.

        Note: if ``use_batch_norm`` is also specified, both will be applied (first
        :class:`~torch.nn.BatchNorm1d`, then :class:`~torch.nn.LayerNorm`).
    use_size_factor_key
        If ``True``, use the :attr:`~anndata.AnnData.obs` column as defined by the
        ``size_factor_key`` parameter in the model's ``setup_anndata`` method as the scaling
        factor in the mean of the conditional distribution. Takes priority over
        ``use_observed_lib_size``.
    use_observed_lib_size
        If ``True``, use the observed library size for RNA as the scaling factor in the mean of the
        conditional distribution.
    library_log_means
        :class:`~numpy.ndarray` of shape ``(1, n_batch)`` of means of the log library sizes that
        parameterize the prior on library size if ``use_size_factor_key`` is ``False`` and
        ``use_observed_lib_size`` is ``False``.
    library_log_vars
        :class:`~numpy.ndarray` of shape ``(1, n_batch)`` of variances of the log library sizes
        that parameterize the prior on library size if ``use_size_factor_key`` is ``False`` and
        ``use_observed_lib_size`` is ``False``.
    extra_decoder_kwargs
        Additional keyword arguments passed into :class:`~scvi.nn.DecoderSCVI`.
    batch_embedding_kwargs
        Keyword arguments passed into :class:`~scvi.nn.Embedding` if ``batch_representation`` is
        set to ``"embedding"``.

    Notes
    -----
    Lifecycle: argument ``batch_representation`` is experimental in v1.2.
    """

    def __init__(
        self,
        n_input: int,
        ##############################
        n_output_niche: int,
        ##############################
        n_batch: int = 0,
        n_labels: int = 0,
        n_hidden: int = 128,
        n_latent: int = 10,
        n_layers: int = 1,
        ##############################
        n_layers_niche: int = 1,
        n_layers_compo: int = 1,
        n_hidden_niche: int = 128,
        n_hidden_compo: int = 128,
        ##############################
        n_continuous_cov: int = 0,
        n_cats_per_cov: list[int] | None = None,
        dropout_rate: float = 0.1,
        dispersion: Literal["gene", "gene-batch", "gene-label", "gene-cell"] = "gene",
        log_variational: bool = True,
        gene_likelihood: Literal["zinb", "nb", "poisson"] = "poisson",
        latent_distribution: Literal["normal", "ln"] = "normal",
        ##############################
        niche_likelihood: Literal["poisson", "gaussian"] = "gaussian",
        cell_rec_weight: float = 1.0,
        latent_kl_weight: float = 1.0,
        spatial_weight: float = 10,
        ##############################
        prior_mixture: bool = False,
        prior_mixture_k: int = 20,
        semisupervised: bool = True,
        linear_classifier: bool = True,
        ##############################
        inpute_covariates_niche_decoder: bool = True,
        encode_covariates: bool = False,
        deeply_inject_covariates: bool = True,
        batch_representation: Literal["one-hot", "embedding"] = "one-hot",
        use_batch_norm: Literal["encoder", "decoder", "none", "both"] = "none",
        use_layer_norm: Literal["encoder", "decoder", "none", "both"] = "both",
        use_size_factor_key: bool = False,
        use_observed_lib_size: bool = True,
        library_log_means: np.ndarray | None = None,
        library_log_vars: np.ndarray | None = None,
        batch_embedding_kwargs: dict | None = None,
        extra_decoder_kwargs: dict | None = None,
        extra_encoder_kwargs: dict | None = None,
        **vae_kwargs,
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
            batch_representation=batch_representation,
            use_size_factor_key=use_size_factor_key,
            use_observed_lib_size=use_observed_lib_size,
            library_log_means=library_log_means,
            library_log_vars=library_log_vars,
            extra_decoder_kwargs=extra_decoder_kwargs,
            batch_embedding_kwargs=batch_embedding_kwargs,
            extra_encoder_kwargs=extra_encoder_kwargs,
            use_batch_norm=use_batch_norm,
            use_layer_norm=use_layer_norm,
            **vae_kwargs,
        )

        self.latent_kl_weight = latent_kl_weight
        self.cell_rec_weight = cell_rec_weight
        self.spatial_weight = spatial_weight
        self.n_output_niche = n_output_niche
        self.niche_likelihood = niche_likelihood
        self.prior_mixture = prior_mixture
        self.semisupervised = semisupervised

        self.batch_representation = batch_representation
        if self.batch_representation == "embedding":
            self.init_embedding(REGISTRY_KEYS.BATCH_KEY, n_batch, **(batch_embedding_kwargs or {}))
            batch_dim = self.get_embedding(REGISTRY_KEYS.BATCH_KEY).embedding_dim
        elif self.batch_representation != "one-hot":
            raise ValueError("`batch_representation` must be one of 'one-hot', 'embedding'.")

        use_batch_norm_encoder = use_batch_norm == "encoder" or use_batch_norm == "both"
        use_batch_norm_decoder = use_batch_norm == "decoder" or use_batch_norm == "both"
        use_layer_norm_encoder = use_layer_norm == "encoder" or use_layer_norm == "both"
        use_layer_norm_decoder = use_layer_norm == "decoder" or use_layer_norm == "both"

        if self.prior_mixture is True:
            if self.semisupervised:
                prior_mixture_k = n_labels
                self.prior_mixture_k = prior_mixture_k

                self.prior_means = torch.nn.Parameter(torch.zeros([prior_mixture_k, n_latent]))
                self.prior_log_scales = torch.nn.Parameter(
                    torch.zeros([prior_mixture_k, n_latent])
                )
                self.prior_logits = torch.nn.Parameter(torch.ones([prior_mixture_k]))

            else:
                self.prior_mixture_k = prior_mixture_k
                self.prior_means = torch.nn.Parameter(torch.randn([prior_mixture_k, n_latent]))
                self.prior_log_scales = torch.nn.Parameter(
                    torch.zeros([prior_mixture_k, n_latent]) - 1.0
                )
                self.prior_logits = torch.nn.Parameter(torch.ones([prior_mixture_k]))

        n_input_encoder = n_input + n_continuous_cov * encode_covariates
        if self.batch_representation == "embedding":
            n_input_encoder += batch_dim * encode_covariates
            cat_list = list([] if n_cats_per_cov is None else n_cats_per_cov)
        else:
            cat_list = [n_batch] + list([] if n_cats_per_cov is None else n_cats_per_cov)

        encoder_cat_list = cat_list if encode_covariates else None
        _extra_encoder_kwargs = extra_encoder_kwargs or {}
        self.z_encoder = Encoder(
            n_input_encoder,
            n_latent,
            n_cat_list=encoder_cat_list,
            n_layers=n_layers,
            n_hidden=n_hidden,
            dropout_rate=dropout_rate,
            distribution=latent_distribution,
            inject_covariates=deeply_inject_covariates,
            use_batch_norm=use_batch_norm_encoder,
            use_layer_norm=use_layer_norm_encoder,
            return_dist=True,
            **_extra_encoder_kwargs,
        )

        n_input_decoder = n_latent + n_continuous_cov
        if self.batch_representation == "embedding":
            n_input_decoder += batch_dim

        _extra_decoder_kwargs = extra_decoder_kwargs or {}

        self.niche_decoder = NicheDecoder(
            n_input=n_input_decoder,
            n_output=n_output_niche,
            n_niche_components=n_labels,
            n_cat_list=cat_list if inpute_covariates_niche_decoder else None,
            n_layers=n_layers_niche,
            n_hidden=n_hidden_niche,
            inject_covariates=deeply_inject_covariates,
            use_batch_norm=use_batch_norm_decoder,
            use_layer_norm=use_layer_norm_decoder,
            dropout_rate=dropout_rate,
            **_extra_decoder_kwargs,
        )

        self.composition_decoder = DirichletDecoder(
            n_input_decoder,
            n_labels,
            n_cat_list=None,  # do not batch-correct the cell type proportions.
            n_layers=n_layers_compo,
            n_hidden=n_hidden_compo,
            inject_covariates=deeply_inject_covariates,
            use_batch_norm=use_batch_norm_decoder,
            use_layer_norm=use_layer_norm_decoder,
            **_extra_decoder_kwargs,
        )

        if self.semisupervised:
            # Classifier takes n_latent as input
            cls_parameters = {
                "n_layers": 0 if linear_classifier else n_layers,
                "n_hidden": 0 if linear_classifier else n_hidden,
                "dropout_rate": dropout_rate,
                "logits": True,  # no Softmax
            }
            self.classifier = Classifier(
                n_latent,
                n_labels=n_labels,
                use_batch_norm=use_batch_norm_encoder,
                use_layer_norm=use_layer_norm_encoder,
                **cls_parameters,
            )
        else:
            self.classifier = None

    @auto_move_data
    def generative(
        self,
        z: torch.Tensor,
        library: torch.Tensor,
        batch_index: torch.Tensor,
        cont_covs: torch.Tensor | None = None,
        cat_covs: torch.Tensor | None = None,
        size_factor: torch.Tensor | None = None,
        y: torch.Tensor | None = None,
        transform_batch: torch.Tensor | None = None,
    ) -> dict[str, Distribution | None]:
        """Run the generative process."""
        from torch.distributions import Categorical, Independent, MixtureSameFamily, Normal
        from torch.nn.functional import linear

        from scvi.distributions import (
            NegativeBinomial,
            Poisson,
            ZeroInflatedNegativeBinomial,
        )

        # TODO: refactor forward function to not rely on y
        # Likelihood distribution
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
            categorical_input = ()

        if transform_batch is not None:
            batch_index = torch.ones_like(batch_index) * transform_batch

        if not self.use_size_factor_key:
            size_factor = library

        if self.batch_representation == "embedding":
            batch_rep = self.compute_embedding(REGISTRY_KEYS.BATCH_KEY, batch_index)
            decoder_input = torch.cat([decoder_input, batch_rep], dim=-1)
            px_scale, px_r, px_rate, px_dropout = self.decoder(
                self.dispersion,
                decoder_input,
                size_factor,
                *categorical_input,
                y,
            )
        else:
            px_scale, px_r, px_rate, px_dropout = self.decoder(
                self.dispersion,
                decoder_input,
                size_factor,
                batch_index,
                *categorical_input,
                y,
            )

        if self.dispersion == "gene-label":
            px_r = linear(
                one_hot(y, self.n_labels), self.px_r
            )  # px_r gets transposed - last dimension is nb genes
        elif self.dispersion == "gene-batch":
            px_r = linear(one_hot(batch_index, self.n_batch), self.px_r)
        elif self.dispersion == "gene":
            px_r = self.px_r

        px_r = torch.exp(px_r)

        if self.gene_likelihood == "zinb":
            px = ZeroInflatedNegativeBinomial(
                mu=px_rate,
                theta=px_r,
                zi_logits=px_dropout,
                scale=px_scale,
            )
        elif self.gene_likelihood == "nb":
            px = NegativeBinomial(mu=px_rate, theta=px_r, scale=px_scale)
        elif self.gene_likelihood == "poisson":
            px = Poisson(px_rate, scale=px_scale)

        # Priors
        if self.use_observed_lib_size:
            pl = None
        else:
            (
                local_library_log_means,
                local_library_log_vars,
            ) = self._compute_local_library_params(batch_index)
            pl = Normal(local_library_log_means, local_library_log_vars.sqrt())

        if self.prior_mixture is True:
            u_prior_logits = self.prior_logits
            u_prior_means = self.prior_means
            u_prior_scales = torch.exp(self.prior_log_scales) + 1e-4

            if self.semisupervised:
                logits_input = (
                    torch.stack(
                        [
                            torch.nn.functional.one_hot(y_i, self.n_labels)
                            if y_i < self.n_labels
                            else torch.zeros(self.n_labels)
                            for y_i in y.ravel()
                        ]
                    )
                    .to(z.device)
                    .float()
                )
                u_prior_logits = u_prior_logits + 10 * logits_input
                u_prior_means = u_prior_means.expand(y.shape[0], -1, -1)
                u_prior_scales = u_prior_scales.expand(y.shape[0], -1, -1)
            cats = Categorical(logits=u_prior_logits)
            normal_dists = Independent(
                Normal(u_prior_means, u_prior_scales), reinterpreted_batch_ndims=1
            )
            pz = MixtureSameFamily(cats, normal_dists)
        else:
            pz = Normal(torch.zeros_like(z), torch.ones_like(z))

        niche_composition = self.composition_decoder(
            decoder_input, batch_index, *categorical_input
        )  # DirichletDecoder, niche_composition is a distribution

        niche_mean, niche_variance = self.niche_decoder(
            decoder_input, batch_index, *categorical_input
        )

        if self.niche_likelihood == "poisson":
            niche_expression = torch.distributions.Poisson(niche_variance)

        else:
            niche_expression = Normal(niche_mean, niche_variance)

        return {
            MODULE_KEYS.PX_KEY: px,
            MODULE_KEYS.PL_KEY: pl,
            MODULE_KEYS.PZ_KEY: pz,
            SCVIVA_MODULE_KEYS.NICHE_MEAN: niche_mean,
            SCVIVA_MODULE_KEYS.NICHE_VARIANCE: niche_variance,
            SCVIVA_MODULE_KEYS.P_NICHE_EXPRESSION: niche_expression,
            SCVIVA_MODULE_KEYS.P_NICHE_COMPOSITION: niche_composition,
        }

    def loss(
        self,
        tensors: dict[str, torch.Tensor],
        inference_outputs: dict[str, torch.Tensor | Distribution | None],
        generative_outputs: dict[str, torch.Tensor | Distribution | None],
        kl_weight: float = 1.0,
        classification_ratio=50,
        epsilon: float = 1e-6,
        n_samples_mixture: int = 10,
    ) -> NicheLossOutput:
        """Compute the loss."""
        from torch.distributions import kl_divergence

        x = tensors[REGISTRY_KEYS.X_KEY]
        if self.semisupervised:
            y = tensors[REGISTRY_KEYS.LABELS_KEY].ravel().long()
            z_mean = inference_outputs[MODULE_KEYS.QZ_KEY].loc
            y_ct = self.classifier(z_mean)
            classification_loss = torch.nn.functional.cross_entropy(y_ct, y, reduction="none")

        if self.prior_mixture is True:
            z = inference_outputs[MODULE_KEYS.QZ_KEY].rsample(
                sample_shape=(n_samples_mixture,)
            )  # sample multiple times
            # sample x n_obs x n_latent
            kl_divergence_z = (
                inference_outputs[MODULE_KEYS.QZ_KEY].log_prob(z).sum(-1)
                - generative_outputs[MODULE_KEYS.PZ_KEY].log_prob(z)
            ).mean(0)

        else:
            kl_divergence_z = kl_divergence(
                inference_outputs[MODULE_KEYS.QZ_KEY],
                generative_outputs[MODULE_KEYS.PZ_KEY],
            ).sum(dim=-1)

        if not self.use_observed_lib_size:
            kl_divergence_l = kl_divergence(
                inference_outputs[MODULE_KEYS.QL_KEY],
                generative_outputs[MODULE_KEYS.PL_KEY],
            ).sum(dim=1)
        else:
            kl_divergence_l = torch.zeros_like(kl_divergence_z)

        reconst_loss_cell = -generative_outputs[MODULE_KEYS.PX_KEY].log_prob(x).sum(-1)

        if self.semisupervised:
            reconst_loss_cell = reconst_loss_cell + classification_ratio * classification_loss

        kl_local_for_warmup = kl_divergence_z
        kl_local_no_warmup = kl_divergence_l

        weighted_kl_local = kl_weight * kl_local_for_warmup + kl_local_no_warmup

        niche_weights = tensors[SCVIVA_REGISTRY_KEYS.NICHE_COMPOSITION_KEY]
        niche_weights = (niche_weights > 0).float()

        z1_mean_niche = tensors[
            SCVIVA_REGISTRY_KEYS.Z1_MEAN_CT_KEY
        ]  # batch times cell_types times n_latent

        reconst_loss_niche = (
            -generative_outputs[SCVIVA_MODULE_KEYS.P_NICHE_EXPRESSION]
            .log_prob(z1_mean_niche)
            .sum(dim=(-1))
        )

        masked_reconst_loss_niche = (reconst_loss_niche * niche_weights).sum(dim=-1)

        true_niche_composition = tensors[SCVIVA_REGISTRY_KEYS.NICHE_COMPOSITION_KEY] + epsilon
        true_niche_composition = true_niche_composition / true_niche_composition.sum(
            dim=-1,
            keepdim=True,
        )

        reconst_niche_composition = generative_outputs[SCVIVA_MODULE_KEYS.P_NICHE_COMPOSITION]

        composition_loss = -reconst_niche_composition.log_prob(true_niche_composition)

        _weighted_reconst_loss_cell = self.cell_rec_weight * reconst_loss_cell
        _weighted_reconst_loss_niche = self.spatial_weight * masked_reconst_loss_niche
        _weighted_composition_loss = self.spatial_weight * composition_loss
        _weighted_kl_local = self.latent_kl_weight * weighted_kl_local

        loss = torch.mean(
            _weighted_reconst_loss_cell
            + _weighted_reconst_loss_niche
            + _weighted_kl_local
            + _weighted_composition_loss
        )

        return NicheLossOutput(
            loss=loss,
            reconstruction_loss=reconst_loss_cell,
            classification_loss=classification_loss.mean() if self.semisupervised else None,
            true_labels=y if self.semisupervised else None,
            logits=y_ct if self.semisupervised else None,
            kl_local={
                MODULE_KEYS.KL_L_KEY: kl_divergence_l,
                MODULE_KEYS.KL_Z_KEY: kl_divergence_z,
            },
            composition_loss=composition_loss,
            niche_loss=masked_reconst_loss_niche,
            extra_metrics={
                SCVIVA_MODULE_KEYS.NLL_NICHE_COMPOSITION_KEY: torch.mean(composition_loss),
                SCVIVA_MODULE_KEYS.NLL_NICHE_EXPRESSION_KEY: torch.mean(masked_reconst_loss_niche),
            },
        )


@dataclass
class NicheLossOutput(LossOutput):
    """Modify loss output to record niche losses."""

    composition_loss: LossRecord | None = None
    niche_loss: LossRecord | None = None

    def __post_init__(self):
        super().__post_init__()

        default = 0 * self.loss
        if self.composition_loss is None:
            object.__setattr__(self, "composition_loss", default)
        if self.niche_loss is None:
            object.__setattr__(self, "niche_loss", default)

        object.__setattr__(self, "composition_loss", self._as_dict("composition_loss"))
        object.__setattr__(self, "niche_loss", self._as_dict("niche_loss"))
