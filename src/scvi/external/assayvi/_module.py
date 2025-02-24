from __future__ import annotations

import logging
from typing import TYPE_CHECKING

import torch
from torch import distributions
from torch.nn.functional import one_hot

from scvi import REGISTRY_KEYS
from scvi.data._constants import ADATA_MINIFY_TYPE
from scvi.module._classifier import Classifier
from scvi.module._constants import MODULE_KEYS
from scvi.module.base import (
    BaseMinifiedModeModuleClass,
    EmbeddingModuleMixin,
    LossOutput,
    MogPrior,
    StandardPrior,
    VampPrior,
    auto_move_data,
)
from scvi.utils import unsupported_if_adata_minified

if TYPE_CHECKING:
    from collections.abc import Callable
    from typing import Literal

    from torch.distributions import Distribution

logger = logging.getLogger(__name__)


class ASSAYVAE(EmbeddingModuleMixin, BaseMinifiedModeModuleClass):
    """Variational auto-encoder :cite:p:`Lopez18`.

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
        * ``"gene-assay"``: parameter is constant per gene per assay.
        * ``"gene-cell"``: parameter is constant per gene per cell.
    log_variational
        If ``True``, use :func:`~torch.log1p` on input data before encoding for numerical stability
        (not normalization).
    gene_likelihood
        Distribution to use for reconstruction in the generative process. One of the following:

        * ``"nb"``: :class:`~scvi.distributions.NegativeBinomial`.
        * ``"zinb"``: :class:`~scvi.distributions.ZeroInflatedNegativeBinomial`.
        * ``"poisson"``: :class:`~scvi.distributions.Poisson`.
        * ``"normal"``: :class:`~torch.distributions.Normal`.
    latent_distribution
        Distribution to use for the latent space. One of the following:

        * ``"normal"``: isotropic normal.
        * ``"ln"``: logistic normal with normal params N(0, 1).
    encode_covariates
        If ``True``, covariates are concatenated to gene expression prior to passing through
        the encoder(s). Else, only gene expression is used.
    deeply_inject_covariates
        If ``True`` and ``n_layers > 1``, covariates are concatenated to the outputs of hidden
        layers in the encoder(s) (if ``encoder_covariates`` is ``True``) and the decoder prior to
        passing through the next layer.
    batch_representation
        Method for encoding batch information. One of the following:

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
    var_activation
        Callable used to ensure positivity of the variance of the variational distribution. Passed
        into :class:`~scvi.nn.Encoder`. Defaults to :func:`~torch.exp`.
    extra_encoder_kwargs
        Additional keyword arguments passed into :class:`~scvi.nn.Encoder`.
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
        n_batch: int = 0,
        n_assay: int = 0,
        n_labels: int = 0,
        n_hidden: int = 128,
        n_latent: int = 10,
        n_layers: int = 1,
        n_continuous_cov: int = 0,
        n_cats_per_cov: list[int] | None = None,
        dropout_rate: float = 0.1,
        dispersion: Literal["gene", "gene-batch", "gene-assay", "gene-cell"] = "gene",
        log_variational: bool = True,
        gene_likelihood: Literal["zinb", "nb", "poisson"] = "nb",
        latent_distribution: Literal["normal", "ln"] = "normal",
        encode_covariates: bool = False,
        encode_assay: bool = True,
        deeply_inject_covariates: bool = False,
        batch_representation: Literal["one-hot", "embedding"] = "one-hot",
        use_batch_norm: Literal["encoder", "decoder", "none", "both"] = "none",
        use_layer_norm: Literal["encoder", "decoder", "none", "both"] = "both",
        var_activation: Callable[[torch.Tensor], torch.Tensor] = None,
        extra_encoder_kwargs: dict | None = None,
        extra_decoder_kwargs: dict | None = None,
        batch_embedding_kwargs: dict | None = None,
        conditional_norm: bool = True,
        conditional_output: bool = True,
        prior: str | None = None,
        pseudoinput_data: dict | None = None,
        n_prior_components: int | None = 30,
        mmd_kernel: str = "rbf",
    ):
        from scvi.nn import DecoderSCVI, Encoder

        super().__init__()

        self.dispersion = dispersion
        self.n_latent = n_latent
        self.log_variational = log_variational
        self.gene_likelihood = gene_likelihood
        self.n_batch = n_batch
        self.n_assay = n_assay
        self.n_labels = n_labels
        self.latent_distribution = latent_distribution
        self.encode_covariates = encode_covariates
        self.use_observed_lib_size = True


        if self.dispersion == "gene":
            self.px_r = torch.nn.Parameter(3.*torch.ones(n_input))
        elif self.dispersion == "gene-batch":
            self.px_r = torch.nn.Parameter(3.*torch.ones(n_input, n_batch))
        elif self.dispersion == "gene-assay":
            self.px_r = torch.nn.Parameter(3.*torch.ones(n_input, n_assay))
        elif self.dispersion == "gene-cell":
            pass
        else:
            raise ValueError(
                "`dispersion` must be one of 'gene', 'gene-batch', 'gene-assay', 'gene-cell'."
            )

        self.batch_representation = batch_representation
        n_cats_per_cov_ = list([] if n_cats_per_cov is None else n_cats_per_cov)
        n_continuous = n_continuous_cov

        if self.batch_representation == "embedding":
            self.init_embedding(REGISTRY_KEYS.BATCH_KEY, n_batch, **(batch_embedding_kwargs or {}))
            n_continuous += self.get_embedding_dim(REGISTRY_KEYS.BATCH_KEY)
        elif self.batch_representation != "one-hot":
            raise ValueError("`batch_representation` must be one of 'one-hot', 'embedding'.")

        use_batch_norm_encoder = use_batch_norm == "encoder" or use_batch_norm == "both"
        use_batch_norm_decoder = use_batch_norm == "decoder" or use_batch_norm == "both"
        use_layer_norm_encoder = use_layer_norm == "encoder" or use_layer_norm == "both"
        use_layer_norm_decoder = use_layer_norm == "decoder" or use_layer_norm == "both"

        if self.batch_representation == "embedding":
            cat_list = n_cats_per_cov_
        else:
            cat_list = [n_batch] + n_cats_per_cov_
        self.encode_assay = encode_assay
        self.batch_representation_encoder = False
        conditional_category = 0
        if encode_assay:
            encode_assay_list = [n_assay]
        else:
            encode_assay_list = [0]
        if not encode_covariates:
            encoder_cat_list = encode_assay_list
            n_cont_encoder = 0
        else:
            if not conditional_norm and self.batch_representation == "embedding":
                # If we are not using conditional norm, we need to pass the batch information
                self.batch_representation_encoder = True
                encoder_cat_list = encode_assay_list + n_cats_per_cov_
                n_cont_encoder = n_continuous
            else:
                encoder_cat_list = encode_assay_list + [n_batch] + n_cats_per_cov_
                n_cont_encoder = n_continuous_cov
            if conditional_norm and not encode_assay:
                conditional_category = 1
        _extra_encoder_kwargs = extra_encoder_kwargs or {}
        self.z_encoder = Encoder(
            n_input,
            n_latent,
            n_continuous=n_cont_encoder,
            n_cat_list=encoder_cat_list,
            n_layers=n_layers,
            n_hidden=n_hidden,
            dropout_rate=dropout_rate,
            distribution=latent_distribution,
            inject_covariates=deeply_inject_covariates,
            use_batch_norm=use_batch_norm_encoder,
            use_layer_norm=use_layer_norm_encoder,
            var_activation=var_activation,
            return_dist=True,
            conditional_norm=conditional_norm,
            conditional_category=conditional_category,
            **_extra_encoder_kwargs,
        )

        _extra_decoder_kwargs = extra_decoder_kwargs or {}
        self.decoder = DecoderSCVI(
            n_latent,
            n_input,
            n_cat_list=cat_list,
            n_continuous=n_continuous,
            n_layers=n_layers,
            n_hidden=n_hidden,
            n_conditions_output=self.n_assay if conditional_output else 0,
            inject_covariates=deeply_inject_covariates,
            use_batch_norm=use_batch_norm_decoder,
            use_layer_norm=use_layer_norm_decoder,
            scale_activation="softmax",
            **_extra_decoder_kwargs,
        )

        cls_parameters = {
            "n_layers": 1,
            "n_hidden": 128,
            "dropout_rate": 0.0,
            "logits": True,
        }

        if self.n_labels > 1:
            self.classifier = Classifier(
                n_latent,
                n_labels=self.n_labels,
                use_batch_norm=False,
                use_layer_norm=True,
                **cls_parameters,
            )
        if prior == "gaussian":
            self.prior = StandardPrior()
        elif prior == "vamp":
            assert pseudoinput_data is not None, (
                "Pseudoinput data must be specified if using VampPrior"
            )
            pseudoinput_data = self._get_inference_input(
                pseudoinput_data,
                full_forward_pass=True
            )
            print('include training.')
            cat_list = [n_batch] + n_cats_per_cov_ + encode_assay_list
            self.prior = VampPrior(
                n_components=n_prior_components,
                inference=self._regular_inference,
                encoder=self.z_encoder,
                pseudoinputs=pseudoinput_data,
                n_cat_list=cat_list,
                trainable_priors=True,
                additional_categorical_covariates=["assay_index"]
            )
        elif prior == "mog":
            self.prior = MogPrior(
                n_components=n_prior_components,
                n_latent=n_latent,
            )
        elif prior == "mog_celltype":
            self.prior = MogPrior(
                n_components=n_labels,
                n_latent=n_latent,
                celltype_bias=True
            )
        else:
            raise ValueError(
                "`prior` must be one of 'gaussian', 'vamp', 'mog', 'mog_celltype'.")

    def _get_inference_input(
        self,
        tensors: dict[str, torch.Tensor | None],
        full_forward_pass: bool = False,
    ) -> dict[str, torch.Tensor | None]:
        """Get input tensors for the inference process."""
        if full_forward_pass or self.minified_data_type is None:
            loader = "full_data"
        elif self.minified_data_type in [
            ADATA_MINIFY_TYPE.LATENT_POSTERIOR,
            ADATA_MINIFY_TYPE.LATENT_POSTERIOR_WITH_COUNTS,
        ]:
            loader = "minified_data"
        else:
            raise NotImplementedError(f"Unknown minified-data type: {self.minified_data_type}")

        if loader == "full_data":
            return {
                MODULE_KEYS.X_KEY: tensors[REGISTRY_KEYS.X_KEY],
                MODULE_KEYS.BATCH_INDEX_KEY: tensors[REGISTRY_KEYS.BATCH_KEY],
                MODULE_KEYS.ASSAY_INDEX_KEY: tensors.get(REGISTRY_KEYS.ASSAY_KEY, None),
                MODULE_KEYS.CONT_COVS_KEY: tensors.get(REGISTRY_KEYS.CONT_COVS_KEY, None),
                MODULE_KEYS.CAT_COVS_KEY: tensors.get(REGISTRY_KEYS.CAT_COVS_KEY, None),
            }
        else:
            return {
                MODULE_KEYS.QZM_KEY: tensors[REGISTRY_KEYS.LATENT_QZM_KEY],
                MODULE_KEYS.QZV_KEY: tensors[REGISTRY_KEYS.LATENT_QZV_KEY],
                REGISTRY_KEYS.OBSERVED_LIB_SIZE: tensors[REGISTRY_KEYS.OBSERVED_LIB_SIZE],
            }

    def _get_generative_input(
        self,
        tensors: dict[str, torch.Tensor],
        inference_outputs: dict[str, torch.Tensor | Distribution | None],
    ) -> dict[str, torch.Tensor | None]:
        """Get input tensors for the generative process."""
        return {
            MODULE_KEYS.Z_KEY: inference_outputs[MODULE_KEYS.Z_KEY],
            MODULE_KEYS.LIBRARY_KEY: inference_outputs[MODULE_KEYS.LIBRARY_KEY],
            MODULE_KEYS.BATCH_INDEX_KEY: tensors[REGISTRY_KEYS.BATCH_KEY],
            MODULE_KEYS.ASSAY_INDEX_KEY: tensors.get(REGISTRY_KEYS.ASSAY_KEY, None),
            MODULE_KEYS.CONT_COVS_KEY: tensors.get(REGISTRY_KEYS.CONT_COVS_KEY, None),
            MODULE_KEYS.CAT_COVS_KEY: tensors.get(REGISTRY_KEYS.CAT_COVS_KEY, None),
        }

    @auto_move_data
    def _regular_inference(
        self,
        x: torch.Tensor,
        batch_index: torch.Tensor,
        assay_index: torch.Tensor | None = None,
        cont_covs: torch.Tensor | None = None,
        cat_covs: torch.Tensor | None = None,
        n_samples: int = 1,
    ) -> dict[str, torch.Tensor | Distribution | None]:
        """Run the regular inference process."""
        x_ = x
        if self.use_observed_lib_size:
            library = torch.log(x.sum(1)).unsqueeze(1)
        if self.log_variational:
            x_ = x_/x_.mean(1).unsqueeze(1)
            x_ = torch.log1p(x_)

        if cat_covs is not None and self.encode_covariates:
            categorical_input = torch.split(cat_covs, 1, dim=1)
        else:
            categorical_input = ()

        if self.encode_covariates and self.batch_representation_encoder:
            batch_rep = self.compute_embedding(REGISTRY_KEYS.BATCH_KEY, batch_index)
            if cont_covs is not None:
                cont_input = torch.cat([cont_covs, batch_rep], dim=-1)
            else:
                cont_input = batch_rep
        else:
            cont_input = cont_covs
        if not self.encode_assay:
            assay_index = None
        else:
            assay_index = assay_index.long()
        qz, z = self.z_encoder(x_, assay_index, batch_index, *categorical_input, cont_input=cont_input)

        if n_samples > 1:
            untran_z = qz.sample((n_samples,))
            z = self.z_encoder.z_transformation(untran_z)
            library = library.unsqueeze(0).expand(
                (n_samples, library.size(0), library.size(1))
            )

        return {
            MODULE_KEYS.Z_KEY: z,
            MODULE_KEYS.QZ_KEY: qz,
            MODULE_KEYS.LIBRARY_KEY: library,
        }

    @auto_move_data
    def generative(
        self,
        z: torch.Tensor,
        library: torch.Tensor,
        batch_index: torch.Tensor,
        assay_index: torch.Tensor | None = None,
        cont_covs: torch.Tensor | None = None,
        cat_covs: torch.Tensor | None = None,
        size_factor: torch.Tensor | None = None, # Consistency
        y: torch.Tensor | None = None,
        transform_batch: torch.Tensor | None = None,
        transform_assay: torch.Tensor | None = None,
    ) -> dict[str, Distribution | None]:
        """Run the generative process."""
        from torch.nn.functional import linear

        from scvi.distributions import (
            NegativeBinomial,
            Normal,
            Poisson,
            ZeroInflatedNegativeBinomial,
        )

        if cat_covs is not None:
            categorical_input = torch.split(cat_covs, 1, dim=1)
        else:
            categorical_input = ()

        if transform_batch is not None:
            batch_index = torch.ones_like(batch_index) * transform_batch
        if transform_assay is not None:
            assay_index = torch.ones_like(assay_index) * transform_assay

        if self.batch_representation == "embedding":
            batch_rep = self.compute_embedding(REGISTRY_KEYS.BATCH_KEY, batch_index)
            if cont_covs is not None:
                cont_input = torch.cat([cont_covs, batch_rep], dim=-1)
            else:
                cont_input = batch_rep
        else:
            cont_input = cont_covs

        px_scale, px_r, px_rate, px_dropout = self.decoder(
            self.dispersion,
            z,
            library,
            batch_index,
            *categorical_input,
            y,
            cont_input=cont_input,
            output_condition=assay_index.long(),
        )

        if self.dispersion == "gene-assay":
            px_r = linear(
                one_hot(assay_index.squeeze(-1).long(), self.n_assay).float(), self.px_r
            )  # px_r gets transposed - last dimension is nb genes
        elif self.dispersion == "gene-batch":
            px_r = linear(one_hot(batch_index.squeeze(-1), self.n_batch).float(), self.px_r)
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
            px = Poisson(rate=px_rate, scale=px_scale)
        elif self.gene_likelihood == "normal":
            px = Normal(px_rate, px_r, normal_mu=px_scale)

        # Priors
        if self.use_observed_lib_size:
            pl = None
        else:
            (
                local_library_log_means,
                local_library_log_vars,
            ) = self._compute_local_library_params(batch_index)
            pl = Normal(local_library_log_means, local_library_log_vars.sqrt())
        pz = Normal(torch.zeros_like(z), torch.ones_like(z))

        return {
            MODULE_KEYS.PX_KEY: px,
            MODULE_KEYS.PL_KEY: pl,
            MODULE_KEYS.PZ_KEY: pz,
        }

    @unsupported_if_adata_minified
    def loss(
        self,
        tensors: dict[str, torch.Tensor],
        inference_outputs: dict[str, torch.Tensor | Distribution | None],
        generative_outputs: dict[str, Distribution | None],
        kl_weight: float = 1.0,
        weight_assay_loss: float = 0.0,
        weight_global: float = 1.0,
        weight_kl_sample: float = 1.0,
        classification_ratio: float = 500.,
    ) -> LossOutput:
        """Compute the loss."""
        from torch.distributions import kl_divergence

        x = tensors[REGISTRY_KEYS.X_KEY]
        y = tensors[REGISTRY_KEYS.LABELS_KEY]

        kl_divergence_z = self.prior.kl(
            qz=inference_outputs[MODULE_KEYS.QZ_KEY],
            z=inference_outputs[MODULE_KEYS.Z_KEY],
            labels=y,
        )

        if self.get_embedding_variational(REGISTRY_KEYS.BATCH_KEY, default_value=False):
            pz_sample = self.compute_embedding(
                REGISTRY_KEYS.BATCH_KEY, tensors[REGISTRY_KEYS.BATCH_KEY], return_dist=True)
            qz_sample = distributions.Normal(
                torch.zeros_like(pz_sample.loc), torch.ones_like(pz_sample.scale))
            kl_divergence_sample = kl_divergence(qz_sample, pz_sample).sum(dim=1)
        else:
            kl_divergence_sample = torch.zeros_like(kl_divergence_z)

        assay_index = tensors.get(REGISTRY_KEYS.ASSAY_KEY, None)
        if weight_assay_loss > 0.0 and assay_index is not None:
            assay_loss = self._compute_assay_penalty(
                inference_outputs[MODULE_KEYS.QZ_KEY].loc,
                assay_index
            )
        else:
            assay_loss = 0.0

        reconst_loss = -generative_outputs[MODULE_KEYS.PX_KEY].log_prob(x).sum(-1)

        kl_global = torch.zeros_like(kl_divergence_z)

        if self.gene_likelihood == "zinb":
            zi = generative_outputs[MODULE_KEYS.PX_KEY].zi_logits
            kl_global -= distributions.Exponential(
                10.*torch.ones_like(zi)).log_prob(torch.exp(zi)).sum(-1)
        if self.gene_likelihood == "zinb" or self.gene_likelihood == "nb":
            theta = generative_outputs[MODULE_KEYS.PX_KEY].theta
            kl_global -= distributions.Exponential(
                torch.ones_like(theta)).log_prob(1/theta).sum(-1)

        weighted_kl_local = kl_weight * (kl_divergence_z + weight_kl_sample * kl_divergence_sample)

        loss = torch.mean(
            reconst_loss + weighted_kl_local + weight_assay_loss * assay_loss +
            weight_global * kl_global)

        if self.n_labels > 1:
            logits = self.classifier(inference_outputs[MODULE_KEYS.Z_KEY])
            classification_loss_ = torch.nn.functional.cross_entropy(logits, y.ravel(), reduction="none")
            mask = (y != self.n_labels)
            classification_loss = classification_ratio * torch.masked_select(
                classification_loss_, mask).mean(0)
            loss += torch.mean(classification_loss)

            return LossOutput(
                loss=loss, reconstruction_loss=reconst_loss, kl_local=kl_divergence_z,
                classification_loss=classification_loss, logits=logits, true_labels=y
            )

        return LossOutput(
            loss=loss,
            reconstruction_loss=reconst_loss,
            kl_local={
                MODULE_KEYS.KL_Z_KEY: kl_divergence_z,
                MODULE_KEYS.KL_SAMPLE_KEY: kl_divergence_sample,
            },
        )

    @torch.inference_mode()
    def sample(
        self,
        tensors: dict[str, torch.Tensor],
        n_samples: int = 1,
        max_poisson_rate: float = 1e8,
    ) -> torch.Tensor:
        r"""Generate predictive samples from the posterior predictive distribution.

        The posterior predictive distribution is denoted as :math:`p(\hat{x} \mid x)`, where
        :math:`x` is the input data and :math:`\hat{x}` is the sampled data.

        We sample from this distribution by first sampling ``n_samples`` times from the posterior
        distribution :math:`q(z \mid x)` for a given observation, and then sampling from the
        likelihood :math:`p(\hat{x} \mid z)` for each of these.

        Parameters
        ----------
        tensors
            Dictionary of tensors passed into :meth:`~scvi.module.VAE.forward`.
        n_samples
            Number of Monte Carlo samples to draw from the distribution for each observation.
        max_poisson_rate
            The maximum value to which to clip the ``rate`` parameter of
            :class:`~scvi.distributions.Poisson`. Avoids numerical sampling issues when the
            parameter is very large due to the variance of the distribution.

        Returns
        -------
        Tensor on CPU with shape ``(n_obs, n_vars)`` if ``n_samples == 1``, else
        ``(n_obs, n_vars,)``.
        """
        from scvi.distributions import Poisson

        inference_kwargs = {"n_samples": n_samples}
        _, generative_outputs = self.forward(
            tensors, inference_kwargs=inference_kwargs, compute_loss=False
        )

        dist = generative_outputs[MODULE_KEYS.PX_KEY]
        if self.gene_likelihood == "poisson":
            dist = Poisson(torch.clamp(dist.rate, max=max_poisson_rate))

        # (n_obs, n_vars) if n_samples == 1, else (n_samples, n_obs, n_vars)
        samples = dist.sample()
        # (n_samples, n_obs, n_vars) -> (n_obs, n_vars, n_samples)
        samples = torch.permute(samples, (1, 2, 0)) if n_samples > 1 else samples

        return samples.cpu()

    def _compute_assay_penalty(
            self, params, assay):
        assay = assay.squeeze(-1).long()
        unique = torch.unique(assay)
        pair_penalty = torch.tensor(0., device=assay.device)
        if len(unique) > 1:
            for i in unique:
                pp = self.mmd(params, mask=(assay == i))
                pair_penalty += pp

        return pair_penalty

    def mmd(self, params, mask=None):
        if mask is not None:
            mod_1 = params[mask]
            mod_2 = params[~mask]
        return rbf_kernel(mod_1, mod_2)


@auto_move_data
def rbf_kernel(x, y, gammas=None):
    """
    Compute the RBF kernel between two tensors.

    Parameters
    ----------
    x : torch.Tensor
        Input tensor of shape (N, D).
    y : torch.Tensor
        Input tensor of shape (M, D).
    gammas : list of float or None
        List of gamma values to compute the kernel with.

    Returns
    -------
    kernel_sum : torch.Tensor
        Tensor of shape (N, M), the sum of RBF kernels over all gamma values.
    """
    if gammas is None:
        gammas = [
            1e-10,
            1e-8,
            1e-6,
            1e-4,
            1e-3,
            1e-2,
            1e-1,
            1,
            2,
            5,
            10,
        ]
    kxy = torch.cdist(x, y).pow(2)
    kxx = torch.cdist(x, x).pow(2)
    kyy = torch.cdist(y, y).pow(2)

    kernel_sum_xy = torch.tensor(0.0, device=x.device)
    kernel_sum_xx = torch.tensor(0.0, device=x.device)
    kernel_sum_yy = torch.tensor(0.0, device=x.device)
    for gamma in gammas:
        kernel_sum_xy += torch.exp(-gamma * kxy).mean()
        kernel_sum_xx += torch.exp(-gamma * kxx).mean()
        kernel_sum_yy += torch.exp(-gamma * kyy).mean()

    return (kernel_sum_xx + kernel_sum_yy - 2 * kernel_sum_xy) / len(gammas)
