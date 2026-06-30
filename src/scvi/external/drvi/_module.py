from __future__ import annotations

from typing import TYPE_CHECKING

import torch
from torch import nn
from torch.nn.functional import linear, one_hot

from scvi import REGISTRY_KEYS
from scvi.distributions import (
    NegativeBinomial,
    Normal,
    Poisson,
    ZeroInflatedNegativeBinomial,
)
from scvi.external.drvi._base_components import DecoderDRVI
from scvi.external.drvi._constants import DRVI_MODULE_KEYS
from scvi.external.drvi._distributions import LogNegativeBinomial
from scvi.module import VAE
from scvi.module._constants import MODULE_KEYS
from scvi.module.base import auto_move_data

if TYPE_CHECKING:
    from collections.abc import Iterable
    from typing import Literal

    from torch.distributions import Distribution


def _resolve_activation(activation: str | type[nn.Module]) -> type[nn.Module]:
    """Resolve an activation name (e.g. ``"elu"``) to an ``nn.Module`` subclass.

    A class / callable returning an ``nn.Module`` is returned unchanged, so users can pass a custom
    activation directly.
    """
    ACTIVATIONS: dict[str, type[nn.Module]] = {
        "relu": nn.ReLU,
        "elu": nn.ELU,
    }

    if isinstance(activation, str):
        try:
            return ACTIVATIONS[activation.lower()]
        except KeyError:
            raise ValueError(
                f"Unknown activation '{activation}'. Choose one of {sorted(ACTIVATIONS)} "
                "or pass an nn.Module subclass."
            ) from None
    return activation


def _resolve_mean_activation(
    activation: str | type[nn.Module] | nn.Module | None,
) -> nn.Module:
    """Resolve a ``mean_activation`` spec to the ``nn.Module`` applied to the latent mean ``q_m``.

    ``None`` is the identity (no-op). An ``nn.Module`` subclass / instance is used directly. A
    string names a :mod:`torch.nn` class (e.g. ``"ReLU"``), with an optional positional argument
    after an underscore (e.g. ``"ELU_0.5"`` -> ``nn.ELU(0.5)``).
    """
    if activation is None:
        return nn.Identity()
    if isinstance(activation, nn.Module):
        return activation
    if isinstance(activation, type):
        return activation()

    name, _, arg = activation.removeprefix("nn.").partition("_")
    cls = getattr(nn, name)
    return cls(float(arg)) if arg else cls()


class DRVIModule(VAE):
    """PyTorch module for DRVI: a disentangled-representation VAE with an additive split decoder.

    Subclasses :class:`~scvi.module.VAE`, reusing its encoder, inference, loss (reconstruction +
    ``KL(qz || N(0, I))``) and likelihoods unchanged. The only deviation is the decoder: the latent
    is split into ``n_split_latent`` groups that are decoded independently and aggregated, which is
    what gives DRVI its interpretable latent factors. See :class:`DecoderDRVI`.

    Parameters
    ----------
    n_input
        Number of input genes.
    n_split_latent
        Number of latent splits. ``None`` or ``-1`` splits every latent dimension (``n_latent``),
        which is the default disentangled setting.
    split_method
        Latent-to-split mapping, one of ``"split_mask"`` or ``"split_map"``. See
        :class:`DecoderDRVI`.
    n_split_output
        Per-split projection output width, i.e. the input dimension to each split's decoder body.
        ``"auto"`` (default) uses ``n_latent``. For ``"split_map"`` this is the learned per-split
        projection size; for ``"split_mask"`` it must equal ``n_latent``. See :class:`DecoderDRVI`.
    split_aggregation
        Per-split aggregation, one of ``"mean"`` or ``"logsumexp"``. See :class:`DecoderDRVI`.
    decoder_reuse_weights
        Decoder weight-sharing policy across splits, one of ``"everywhere"`` (default),
        ``"hidden"``, ``"last"``, ``"hidden_except_first"`` or ``"nowhere"``. See
        :class:`DecoderDRVI`.
    n_latent
        Dimensionality of the latent space.
    gene_likelihood
        Reconstruction likelihood. One of:

        * ``"nb"`` (default) - negative binomial (mean ``= library * softmax(decoder)``).
        * ``"pnb"`` - parametrized negative binomial: same mean, but modeled in **log space** via
          :class:`~scvi.external.drvi.LogNegativeBinomial`. This is the numerically stable,
          additive form that composes the per-split log contributions of the decoder
          (recommended with ``split_aggregation="logsumexp"``).
        * ``"zinb"`` - zero-inflated negative binomial.
        * ``"poisson"`` - Poisson.
        * ``"normal"`` - Gaussian with the mean modeled directly (no library/softmax) and per-gene
          variance modeled in log space; use with continuous (e.g. log-normalized) data.
        * ``"normal_unit_var"`` - Gaussian with unit variance.
    activation_fn
        Hidden-layer activation for both the encoder and the split decoder. Either a name (one of
        ``"elu"`` (default, DRVI's choice),``"relu"``, or an ``nn.Module`` subclass.
    mean_activation
        Activation applied to the encoder's latent mean ``q_m``. ``None`` (default) is a no-op; a
        non-identity activation (e.g. ``"ReLU"``) constrains the latent space (e.g. non-negative
        means). Accepts the name of a :mod:`torch.nn` class with an optional positional argument
        after an underscore (e.g. ``"ReLU"``, ``"GELU"``, ``"ELU_0.5"``, ``"LeakyReLU_0.2"``), an
        ``nn.Module`` subclass, or an instance. Applied only to the latent encoder, not the
        library encoder.
    extra_encoder_kwargs
        Extra keyword arguments for the encoder :class:`~scvi.nn.FCLayers`.
    extra_decoder_kwargs
        Extra keyword arguments for the decoder :class:`~scvi.external.drvi.SplitFCLayers`.
    **kwargs
        Additional keyword arguments for :class:`~scvi.module.VAE`.
    """

    def __init__(
        self,
        n_input: int,
        n_split_latent: int | None = None,
        split_method: Literal["split_mask", "split_map"] = "split_map",
        n_split_output: int | Literal["auto"] = "auto",
        split_aggregation: Literal["mean", "logsumexp"] = "logsumexp",
        decoder_reuse_weights: Literal[
            "everywhere", "last", "hidden", "nowhere", "hidden_except_first"
        ] = "everywhere",
        n_latent: int = 128,
        n_hidden: int = 256,
        n_layers: int = 1,
        dropout_rate: float = 0.1,
        n_continuous_cov: int = 0,
        n_cats_per_cov: Iterable[int] | None = None,
        dispersion: Literal["gene", "gene-batch", "gene-label", "gene-cell"] = "gene",
        gene_likelihood: Literal[
            "nb", "pnb", "zinb", "poisson", "normal", "normal_unit_var"
        ] = "pnb",
        deeply_inject_covariates: bool = False,
        use_batch_norm: Literal["encoder", "decoder", "none", "both"] = "none",
        use_layer_norm: Literal["encoder", "decoder", "none", "both"] = "both",
        activation_fn: str | type[nn.Module] = "elu",
        mean_activation: str | type[nn.Module] | nn.Module | None = None,
        extra_encoder_kwargs: dict | None = None,
        extra_decoder_kwargs: dict | None = None,
        **kwargs,
    ):
        # hidden-layer activation for both encoder and (split) decoder; DRVI defaults to ELU.
        activation = _resolve_activation(activation_fn)
        extra_encoder_kwargs = {"activation_fn": activation, **(extra_encoder_kwargs or {})}
        decoder_extra_kwargs = {"activation_fn": activation, **(extra_decoder_kwargs or {})}

        super().__init__(
            n_input,
            n_latent=n_latent,
            n_hidden=n_hidden,
            n_layers=n_layers,
            n_continuous_cov=n_continuous_cov,
            n_cats_per_cov=n_cats_per_cov,
            dispersion=dispersion,
            gene_likelihood=gene_likelihood,
            dropout_rate=dropout_rate,
            deeply_inject_covariates=deeply_inject_covariates,
            use_batch_norm=use_batch_norm,
            use_layer_norm=use_layer_norm,
            extra_encoder_kwargs=extra_encoder_kwargs,
            **kwargs,
        )

        # apply the optional latent-mean activation to the latent (z) encoder.
        mean_act = _resolve_mean_activation(mean_activation)
        if not isinstance(mean_act, nn.Identity):
            self.z_encoder.mean_encoder = nn.Sequential(self.z_encoder.mean_encoder, mean_act)

        if n_split_latent is None or n_split_latent == -1:
            n_split_latent = n_latent
        self.n_split_latent = n_split_latent
        self.split_method = split_method
        self.split_aggregation = split_aggregation
        # analysis flags read by the generative / interpretability mixins
        self.inspect_mode = False
        self.fully_deterministic = False

        # With the embedding batch representation the batch is injected into each split as a
        # continuous covariate (the learned batch embedding) rather than one-hot in n_cat_list.
        if self.batch_representation == "embedding":
            batch_dim = self.get_embedding(REGISTRY_KEYS.BATCH_KEY).embedding_dim
            cat_list = list(n_cats_per_cov or [])
            decoder_n_continuous_cov = n_continuous_cov + batch_dim
        else:
            cat_list = [self.n_batch] + list(n_cats_per_cov or [])
            decoder_n_continuous_cov = n_continuous_cov
        use_batch_norm_decoder = use_batch_norm in ("decoder", "both")
        use_layer_norm_decoder = use_layer_norm in ("decoder", "both")

        # replace the inherited decoder with drvi additive decoder
        self.decoder = DecoderDRVI(
            n_latent,
            n_input,
            n_split=self.n_split_latent,
            split_method=split_method,
            n_split_output=n_split_output,
            split_aggregation=split_aggregation,
            n_cat_list=cat_list,
            n_continuous_cov=decoder_n_continuous_cov,
            n_layers=n_layers,
            n_hidden=n_hidden,
            reuse_weights=decoder_reuse_weights,
            inject_covariates=deeply_inject_covariates,
            use_batch_norm=use_batch_norm_decoder,
            use_layer_norm=use_layer_norm_decoder,
            dropout_rate=dropout_rate,
            model_cell_dispersion=dispersion == "gene-cell",
            model_zero_inflation=gene_likelihood == "zinb",
            **decoder_extra_kwargs,
        )

    def _regular_inference(self, *args, **kwargs):
        # super()'s implementation is @auto_move_data-decorated and handles device placement;
        # we only post-process its output here.
        outputs = super()._regular_inference(*args, **kwargs)
        if self.fully_deterministic:
            qz = outputs[MODULE_KEYS.QZ_KEY]
            outputs[MODULE_KEYS.Z_KEY] = self.z_encoder.z_transformation(qz.loc)
        return outputs

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
    ) -> dict[str, Distribution | torch.Tensor | None]:
        """Run the generative process via the additive decoder.

        Mirrors :meth:`scvi.module.VAE.generative` but passes continuous covariates to the decoder
        separately (so the split transformation only sees the latent dimensions). With the
        embedding batch representation, the learned batch embedding is injected into each split as
        an extra continuous covariate.
        """
        categorical_input = torch.split(cat_covs, 1, dim=1) if cat_covs is not None else ()
        if transform_batch is not None:
            batch_index = torch.ones_like(batch_index) * transform_batch
        if not self.use_size_factor_key:
            size_factor = library
        elif size_factor is None:
            raise ValueError(
                "DRVIModule was initialized with use_size_factor_key=True, but no size_factor "
                "tensor was provided to generative(). Ensure setup_anndata received a "
                "size_factor_key and that the tensor dictionary includes it."
            )

        # batch handling: one-hot batch is injected via the decoder's n_cat_list; the embedding
        # batch representation is concatenated to each split as a continuous covariate.
        if self.batch_representation == "embedding":
            batch_rep = self.compute_embedding(REGISTRY_KEYS.BATCH_KEY, batch_index)
            decoder_cont = (
                batch_rep if cont_covs is None else torch.cat([cont_covs, batch_rep], dim=-1)
            )
            decoder_cats = categorical_input
        else:
            decoder_cont = cont_covs
            decoder_cats = (batch_index, *categorical_input)

        self.decoder.inspect_mode = self.inspect_mode
        # the decoder returns log-space per-gene parameters. Labels are not part of the decoder's
        # n_cat_list, so y is not passed here (unlike scvi's VAE); it is only used for the
        # gene-label dispersion below.
        px_scale_logit, px_r_logit, px_dropout_logit, px_scale_logit_per_split = self.decoder(
            z,
            *decoder_cats,
            cont=decoder_cont,
        )
        # log_softmax over genes, then add the (log) library size to get log(mu)
        px_scale_log = px_scale_logit - torch.logsumexp(px_scale_logit, dim=-1, keepdim=True)
        px_rate_log = size_factor + px_scale_log  # size_factor == log(library size)

        # dispersion logit, before exponentiation: log-theta for the (log-)NB likelihoods and
        # log-variance for the normal likelihoods.
        if self.dispersion == "gene-label":
            px_r_logit = linear(one_hot(y.squeeze(-1), self.n_labels).float(), self.px_r)
        elif self.dispersion == "gene-batch":
            px_r_logit = linear(one_hot(batch_index.squeeze(-1), self.n_batch).float(), self.px_r)
        elif self.dispersion == "gene":
            px_r_logit = self.px_r

        if self.gene_likelihood == "pnb":
            px = LogNegativeBinomial(log_m=px_rate_log, log_r=px_r_logit, log_scale=px_scale_log)
        elif self.gene_likelihood in ("nb", "zinb", "poisson"):
            px_scale = torch.exp(px_scale_log)
            px_rate = torch.exp(px_rate_log)
            if self.gene_likelihood == "nb":
                px = NegativeBinomial(mu=px_rate, theta=torch.exp(px_r_logit), scale=px_scale)
            elif self.gene_likelihood == "zinb":
                px = ZeroInflatedNegativeBinomial(
                    mu=px_rate,
                    theta=torch.exp(px_r_logit),
                    zi_logits=px_dropout_logit,
                    scale=px_scale,
                )
            else:  # poisson
                px = Poisson(rate=px_rate, scale=px_scale)
        elif self.gene_likelihood == "normal":
            # Gaussian with the mean modeled directly (the raw log-space decoder output, no
            # library/softmax) and per-gene variance modeled in log space.
            var = torch.nan_to_num(torch.exp(px_r_logit), posinf=100.0, neginf=0.0) + 1e-8
            px = Normal(px_scale_logit, var.sqrt(), normal_mu=px_scale_logit)
        elif self.gene_likelihood == "normal_unit_var":
            px = Normal(px_scale_logit, torch.ones_like(px_scale_logit), normal_mu=px_scale_logit)
        else:
            raise ValueError(f"Unknown gene_likelihood: {self.gene_likelihood}")

        if self.use_observed_lib_size:
            pl = None
        else:
            local_library_log_means, local_library_log_vars = self._compute_local_library_params(
                batch_index
            )
            pl = Normal(local_library_log_means, local_library_log_vars.sqrt())
        pz = Normal(torch.zeros_like(z), torch.ones_like(z))

        return {
            MODULE_KEYS.PX_KEY: px,
            MODULE_KEYS.PL_KEY: pl,
            MODULE_KEYS.PZ_KEY: pz,
            DRVI_MODULE_KEYS.PX_UNAGGREGATED_PARAMS_KEY: px_scale_logit_per_split,
        }
