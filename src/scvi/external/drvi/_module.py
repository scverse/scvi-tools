from __future__ import annotations

from typing import TYPE_CHECKING

import torch

from scvi.module import VAE
from scvi.module._constants import MODULE_KEYS
from scvi.module.base import auto_move_data

from ._base_components import SplitDecoder
from ._constants import DRVI_MODULE_KEYS
from ._distributions import LogNegativeBinomial

if TYPE_CHECKING:
    from collections.abc import Iterable
    from typing import Literal

    from torch.distributions import Distribution


class DRVIModule(VAE):
    """PyTorch module for DRVI: a disentangled-representation VAE with an additive split decoder.

    Subclasses :class:`~scvi.module.VAE`, reusing its encoder, inference, loss (reconstruction +
    ``KL(qz || N(0, I))``) and likelihoods unchanged. The only deviation is the decoder: the latent
    is split into ``n_split_latent`` groups that are decoded independently and aggregated, which is
    what gives DRVI its interpretable latent factors. See :class:`SplitDecoder`.

    Parameters
    ----------
    n_input
        Number of input genes.
    n_split_latent
        Number of latent splits. ``None`` or ``-1`` splits every latent dimension (``n_latent``),
        which is the default disentangled setting.
    split_method
        Latent-to-split mapping, one of ``"split_diag"`` or ``"split_map"``. See
        :class:`SplitDecoder`.
    split_aggregation
        Per-split aggregation, one of ``"mean"`` or ``"logsumexp"``. See :class:`SplitDecoder`.
    reuse_weights
        Whether decoder weights are shared across splits (``True``) or learned per split
        (``False``, via :class:`stacked_linear.StackedLinearLayer`).
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
    **kwargs
        Additional keyword arguments for :class:`~scvi.module.VAE`.
    """

    def __init__(
        self,
        n_input: int,
        n_split_latent: int | None = None,
        split_method: Literal["split_diag", "split_map"] = "split_map",
        split_aggregation: Literal["mean", "logsumexp"] = "logsumexp",
        reuse_weights: bool = True,
        n_latent: int = 32,
        n_hidden: int = 128,
        n_layers: int = 1,
        n_continuous_cov: int = 0,
        n_cats_per_cov: Iterable[int] | None = None,
        dispersion: Literal["gene", "gene-batch", "gene-label", "gene-cell"] = "gene",
        gene_likelihood: Literal[
            "nb", "pnb", "zinb", "poisson", "normal", "normal_unit_var"
        ] = "nb",
        deeply_inject_covariates: bool = False,
        use_batch_norm: Literal["encoder", "decoder", "none", "both"] = "none",
        use_layer_norm: Literal["encoder", "decoder", "none", "both"] = "both",
        **kwargs,
    ):
        super().__init__(
            n_input,
            n_latent=n_latent,
            n_hidden=n_hidden,
            n_layers=n_layers,
            n_continuous_cov=n_continuous_cov,
            n_cats_per_cov=n_cats_per_cov,
            dispersion=dispersion,
            gene_likelihood=gene_likelihood,
            deeply_inject_covariates=deeply_inject_covariates,
            use_batch_norm=use_batch_norm,
            use_layer_norm=use_layer_norm,
            **kwargs,
        )

        if n_split_latent is None or n_split_latent == -1:
            n_split_latent = n_latent
        self.n_split_latent = n_split_latent
        self.split_method = split_method
        self.split_aggregation = split_aggregation
        # analysis flags read by the generative / interpretability mixins
        self.inspect_mode = False
        self.fully_deterministic = False

        cat_list = [self.n_batch] + list(n_cats_per_cov or [])
        use_batch_norm_decoder = use_batch_norm in ("decoder", "both")
        use_layer_norm_decoder = use_layer_norm in ("decoder", "both")

        # replace the inherited DecoderSCVI with the additive split decoder
        self.decoder = SplitDecoder(
            n_latent,
            n_input,
            n_split=self.n_split_latent,
            split_method=split_method,
            split_aggregation=split_aggregation,
            n_cat_list=cat_list,
            n_continuous_cov=n_continuous_cov,
            n_layers=n_layers,
            n_hidden=n_hidden,
            reuse_weights=reuse_weights,
            inject_covariates=deeply_inject_covariates,
            use_batch_norm=use_batch_norm_decoder,
            use_layer_norm=use_layer_norm_decoder,
            scale_activation="softmax",
            dropout_rate=0.0,
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
    ) -> dict[str, Distribution | None]:
        """Run the generative process via the additive split decoder.

        Mirrors :meth:`scvi.module.VAE.generative` but passes continuous covariates to the decoder
        separately (so the split transformation only sees the latent dimensions) and supports the
        one-hot batch representation only.
        """
        from torch.nn.functional import linear, one_hot

        from scvi.distributions import (
            NegativeBinomial,
            Normal,
            Poisson,
            ZeroInflatedNegativeBinomial,
        )

        categorical_input = torch.split(cat_covs, 1, dim=1) if cat_covs is not None else ()
        if transform_batch is not None:
            batch_index = torch.ones_like(batch_index) * transform_batch
        if not self.use_size_factor_key:
            size_factor = library

        self.decoder.inspect_mode = self.inspect_mode
        px_scale, px_r, px_rate, px_dropout, px_agg = self.decoder(
            self.dispersion,
            z,
            size_factor,
            batch_index,
            *categorical_input,
            y,
            cont=cont_covs,
        )

        # dispersion logit, before exponentiation: log-theta for the (log-)NB likelihoods and
        # log-variance for the normal likelihoods.
        if self.dispersion == "gene-label":
            px_r_logit = linear(one_hot(y.squeeze(-1), self.n_labels).float(), self.px_r)
        elif self.dispersion == "gene-batch":
            px_r_logit = linear(one_hot(batch_index.squeeze(-1), self.n_batch).float(), self.px_r)
        elif self.dispersion == "gene":
            px_r_logit = self.px_r
        else:  # gene-cell: per cell-gene logit produced by the decoder
            px_r_logit = px_r

        if self.gene_likelihood == "zinb":
            px = ZeroInflatedNegativeBinomial(
                mu=px_rate, theta=torch.exp(px_r_logit), zi_logits=px_dropout, scale=px_scale
            )
        elif self.gene_likelihood == "nb":
            px = NegativeBinomial(mu=px_rate, theta=torch.exp(px_r_logit), scale=px_scale)
        elif self.gene_likelihood == "pnb":
            # parametrized negative binomial: model the mean in log space, so the additive
            # (logsumexp) split decoder composes per-split contributions correctly and the NB
            # log-prob stays numerically stable. mu = lib * softmax(agg), theta = exp(px_r_logit).
            # log_softmax over genes
            log_scale = px_agg - torch.logsumexp(px_agg, dim=-1, keepdim=True)
            log_m = library + log_scale  # library == log(observed library size)
            px = LogNegativeBinomial(log_m=log_m, log_r=px_r_logit, log_scale=log_scale)
        elif self.gene_likelihood == "poisson":
            px = Poisson(rate=px_rate, scale=px_scale)
        elif self.gene_likelihood == "normal":
            # Gaussian with the mean modeled directly (no library/softmax) and per-gene variance
            # modeled in log space.
            var = torch.nan_to_num(torch.exp(px_r_logit), posinf=100.0, neginf=0.0) + 1e-8
            px = Normal(px_agg, var.sqrt(), normal_mu=px_agg)
        elif self.gene_likelihood == "normal_unit_var":
            px = Normal(px_agg, torch.ones_like(px_agg), normal_mu=px_agg)
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

        outputs = {
            MODULE_KEYS.PX_KEY: px,
            MODULE_KEYS.PL_KEY: pl,
            MODULE_KEYS.PZ_KEY: pz,
        }
        if self.inspect_mode:
            outputs[DRVI_MODULE_KEYS.PX_UNAGGREGATED_PARAMS_KEY] = {
                "mean": self.decoder._split_scale_cache
            }
        return outputs
