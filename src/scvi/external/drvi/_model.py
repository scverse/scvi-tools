from __future__ import annotations

import logging
from typing import TYPE_CHECKING

from scvi import REGISTRY_KEYS
from scvi.data import AnnDataManager
from scvi.data.fields import (
    CategoricalJointObsField,
    CategoricalObsField,
    LayerField,
    NumericalJointObsField,
)
from scvi.model.base import (
    ArchesMixin,
    BaseModelClass,
    RNASeqMixin,
    UnsupervisedTrainingMixin,
    VAEMixin,
)
from scvi.utils import setup_anndata_dsp

from ._generative_mixin import GenerativeMixin
from ._interpretability_mixin import InterpretabilityMixin
from ._module import DRVIModule

if TYPE_CHECKING:
    from typing import Literal

    from anndata import AnnData

logger = logging.getLogger(__name__)


class DRVI(
    RNASeqMixin,
    VAEMixin,
    ArchesMixin,
    UnsupervisedTrainingMixin,
    BaseModelClass,
    GenerativeMixin,
    InterpretabilityMixin,
):
    """Disentangled Representation Variational Inference :cite:p:`Moinfar2024`.

    DRVI is an unsupervised deep generative model that learns an **interpretable, disentangled**
    latent representation of single-cell omics data. Disentanglement is induced in the decoder: the
    latent is split into independent groups, each decoded separately and aggregated (see
    :class:`~scvi.external.drvi.SplitDecoder`). Built on :class:`~scvi.module.VAE`.

    Parameters
    ----------
    adata
        AnnData object registered via :meth:`~scvi.external.DRVI.setup_anndata`. May be ``None``
        when initializing from a ``registry`` (e.g. out-of-core training with a datamodule).
    registry
        Setup registry (e.g. from a datamodule such as
        :class:`~scvi.dataloaders.AnnbatchDataModule`) to initialize the model without an in-memory
        AnnData. Mutually exclusive with ``adata``.
    n_latent
        Dimensionality of the latent space.
    n_split_latent
        Number of latent splits. ``None`` (default) splits every latent dimension.
    split_method
        Latent-to-split mapping, ``"split_diag"`` or ``"split_map"`` (default).
    split_aggregation
        Per-split aggregation, ``"mean"`` or ``"logsumexp"`` (default).
    **model_kwargs
        Keyword args for :class:`~scvi.external.drvi.DRVIModule`.

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
        split_method: Literal["split_diag", "split_map"] = "split_map",
        split_aggregation: Literal["mean", "logsumexp"] = "logsumexp",
        **model_kwargs,
    ):
        super().__init__(adata, registry)

        if adata is not None:
            n_cats_per_cov = (
                self.adata_manager.get_state_registry(REGISTRY_KEYS.CAT_COVS_KEY).n_cats_per_key
                if REGISTRY_KEYS.CAT_COVS_KEY in self.adata_manager.data_registry
                else None
            )
        else:
            # registry path (e.g. out-of-core training with a datamodule)
            cat_cov_sr = self.registry["field_registries"][REGISTRY_KEYS.CAT_COVS_KEY][
                "state_registry"
            ]
            n_cats_per_cov = tuple(cat_cov_sr["n_cats_per_key"]) if cat_cov_sr else None

        self.module = self._module_cls(
            n_input=self.summary_stats.n_vars,
            n_batch=self.summary_stats.n_batch,
            n_labels=self.summary_stats.get("n_labels", 1),
            n_continuous_cov=self.summary_stats.get("n_extra_continuous_covs", 0),
            n_cats_per_cov=n_cats_per_cov,
            n_latent=n_latent,
            n_split_latent=n_split_latent,
            split_method=split_method,
            split_aggregation=split_aggregation,
            **model_kwargs,
        )

        self._model_summary_string = (
            f"DRVI model with n_latent: {self.module.n_latent}, "
            f"n_split_latent: {self.module.n_split_latent}, "
            f"split_method: '{self.module.split_method}', "
            f"split_aggregation: '{self.module.split_aggregation}'."
        )
        self.init_params_ = self._get_init_params(locals())
        logger.info("The model has been initialized")

    @classmethod
    @setup_anndata_dsp.dedent
    def setup_anndata(
        cls,
        adata: AnnData,
        layer: str | None = None,
        batch_key: str | None = None,
        labels_key: str | None = None,
        categorical_covariate_keys: list[str] | None = None,
        continuous_covariate_keys: list[str] | None = None,
        **kwargs,
    ):
        """Set up :class:`~anndata.AnnData` for DRVI.

        Parameters
        ----------
        %(param_adata)s
        %(param_layer)s
        %(param_batch_key)s
        %(param_labels_key)s
        %(param_cat_cov_keys)s
        %(param_cont_cov_keys)s
        """
        setup_method_args = cls._get_setup_method_args(**locals())
        anndata_fields = [
            LayerField(REGISTRY_KEYS.X_KEY, layer, is_count_data=True),
            CategoricalObsField(REGISTRY_KEYS.BATCH_KEY, batch_key),
            CategoricalObsField(REGISTRY_KEYS.LABELS_KEY, labels_key),
            CategoricalJointObsField(REGISTRY_KEYS.CAT_COVS_KEY, categorical_covariate_keys),
            NumericalJointObsField(REGISTRY_KEYS.CONT_COVS_KEY, continuous_covariate_keys),
        ]
        adata_manager = AnnDataManager(fields=anndata_fields, setup_method_args=setup_method_args)
        adata_manager.register_fields(adata, **kwargs)
        cls.register_manager(adata_manager)
