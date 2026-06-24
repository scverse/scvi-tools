from __future__ import annotations

import logging
from typing import TYPE_CHECKING

from scvi.external.drvi._generative_mixin import GenerativeMixin
from scvi.external.drvi._interpretability_mixin import InterpretabilityMixin
from scvi.external.drvi._module import DRVIModule
from scvi.model import SCVI

if TYPE_CHECKING:
    from typing import Literal

    from anndata import AnnData

logger = logging.getLogger(__name__)


class DRVI(SCVI, GenerativeMixin, InterpretabilityMixin):
    """Disentangled Representation Variational Inference :cite:p:`Moinfar2024`.

    DRVI is an unsupervised deep generative model that learns an **interpretable, disentangled**
    latent representation of single-cell omics data. Disentanglement is induced in the decoder: the
    latent is split into independent groups, each decoded separately and aggregated (see
    :class:`~scvi.external.drvi.DecoderDRVI`).

    DRVI is built directly on :class:`~scvi.model.SCVI`: it reuses SCVI's encoder, likelihoods,
    covariate/batch handling, minified-mode, scArches and out-of-core (datamodule) machinery
    unchanged, and only swaps the decoder for the additive split decoder and adds the
    disentanglement/interpretability logic. It therefore exposes the full SCVI interface
    (:meth:`~scvi.model.SCVI.setup_anndata`, ``minify_adata``, ``get_batch_representation``,
    ``size_factor_key``, …) plus the DRVI-specific knobs below and the interpretability methods
    from :class:`~scvi.external.drvi.InterpretabilityMixin`.

    Parameters
    ----------
    adata
        AnnData object registered via :meth:`~scvi.external.DRVI.setup_anndata`. May be ``None``
        when initializing from a ``registry`` (e.g. out-of-core training with a datamodule) or to
        defer module initialization to ``train`` (datamodule path).
    registry
        Setup registry (e.g. from a datamodule such as
        :class:`~scvi.dataloaders.AnnbatchDataModule`) to initialize without an in-memory AnnData.
    n_latent
        Dimensionality of the latent space.
    n_split_latent
        Number of latent splits. ``None`` (default) splits every latent dimension.
    split_method
        Latent-to-split mapping, ``"split_mask"`` or ``"split_map"`` (default).
    split_aggregation
        Per-split aggregation, ``"mean"`` or ``"logsumexp"`` (default).
    gene_likelihood
        Reconstruction likelihood; ``"pnb"`` (default, log-space parametrized NB), ``"nb"``,
        ``"zinb"``, ``"poisson"``, ``"normal"`` or ``"normal_unit_var"``.
    **kwargs
        Additional keyword args for :class:`~scvi.model.SCVI` /
        :class:`~scvi.external.drvi.DRVIModule`
        (e.g. ``n_hidden``, ``n_layers``, ``n_split_output``, ``dispersion``,
        ``batch_representation``, ``activation_fn``, ``use_observed_lib_size``).

    Examples
    --------
    >>> adata = anndata.read_h5ad(path_to_anndata)
    >>> scvi.external.DRVI.setup_anndata(adata, batch_key="batch")
    >>> model = scvi.external.DRVI(adata, n_latent=32)
    >>> model.train()
    >>> adata.obsm["X_drvi"] = model.get_latent_representation()
    """

    _module_cls = DRVIModule

    def __init__(
        self,
        adata: AnnData | None = None,
        registry: dict | None = None,
        n_latent: int = 32,
        n_split_latent: int | None = None,
        split_method: Literal["split_mask", "split_map"] = "split_map",
        split_aggregation: Literal["mean", "logsumexp"] = "logsumexp",
        gene_likelihood: Literal[
            "nb", "pnb", "zinb", "poisson", "normal", "normal_unit_var"
        ] = "pnb",
        **kwargs,
    ):
        # SCVI.__init__ handles the adata / registry / deferred (datamodule) paths, library size
        # init, size-factor and minified-mode wiring, and builds ``self._module_cls`` (DRVIModule).
        # The DRVI-specific knobs flow through to DRVIModule via ``**kwargs``.
        super().__init__(
            adata,
            registry,
            n_latent=n_latent,
            gene_likelihood=gene_likelihood,
            n_split_latent=n_split_latent,
            split_method=split_method,
            split_aggregation=split_aggregation,
            **kwargs,
        )
        self._model_summary_string = (
            "Model with the following params:\n"
            f"  n_latent: {self.module.n_latent}\n"
            f"  n_split_latent: {self.module.n_split_latent}\n"
            f"  split_method: {self.module.split_method!r}\n"
            f"  split_aggregation: {self.module.split_aggregation!r}"
            f". gene_likelihood: {self.module.gene_likelihood!r}"
        )
        self.init_params_ = self._get_init_params(locals())
        logger.info("The model has been initialized")
