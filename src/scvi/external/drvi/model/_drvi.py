from __future__ import annotations

import logging
from typing import TYPE_CHECKING
from importlib.metadata import version

import numpy as np
from scvi import REGISTRY_KEYS
from scvi.data import AnnDataManager
from scvi.data.fields import CategoricalObsField, LayerField, NumericalJointObsField, CategoricalJointObsField
from scvi.model.base import BaseModelClass, RNASeqMixin, UnsupervisedTrainingMixin, VAEMixin
from scvi.utils import setup_anndata_dsp

from scvi.external.drvi.model.base import DRVIArchesMixin, GenerativeMixin
from scvi.external.drvi.module import DRVIModule

if TYPE_CHECKING:
    from collections.abc import Sequence
    from typing import Any, Literal

    from anndata import AnnData


_DRVI_LATENT_QZM = "_drvi_latent_qzm"
_DRVI_LATENT_QZV = "_drvi_latent_qzv"
_DRVI_OBSERVED_LIB_SIZE = "_drvi_observed_lib_size"

SCVI_VERSION = version("scvi-tools")

logger = logging.getLogger(__name__)


class DRVI(RNASeqMixin, VAEMixin, DRVIArchesMixin, UnsupervisedTrainingMixin, BaseModelClass, GenerativeMixin):
    """DRVI model based on scvi-tools framework for disentangled representation learning.

    Parameters
    ----------
    adata
        AnnData object that has been registered via :meth:`~drvi.model.DRVI.setup_anndata`.
    n_latent
        Dimensionality of the latent space.
    encoder_dims
        Number of nodes in hidden layers of the encoder.
    decoder_dims
        Number of nodes in hidden layers of the decoder.
    prior
        Prior model type.
    categorical_embedding_dims
        Dictionary mapping categorical covariate names to their embedding dimensions.
        Used only if `covariate_modeling_strategy` passed to DRVIModule is based on embedding (not onehot encoding).
        Keys should match the covariate names used in :meth:`~drvi.model.DRVI.setup_anndata`.
        If not provided, default embedding dimension of 10 is used for all covariates.
    **model_kwargs
        Additional keyword arguments passed to :class:`~drvi.model.DRVIModule`.

    Examples
    --------
    >>> adata = anndata.read_h5ad(path_to_anndata)
    >>> drvi.model.DRVI.setup_anndata(adata, categorical_covariate_keys=["batch"])
    >>> vae = drvi.model.DRVI(adata)
    >>> vae.train()
    >>> adata.obsm["latent"] = vae.get_latent_representation()
    """

    _module_cls = DRVIModule
    _LATENT_QZM_KEY = _DRVI_LATENT_QZM
    _LATENT_QZV_KEY = _DRVI_LATENT_QZV
    _DEFAULT_CATEGORICAL_EMBEDDING_DIM = 10

    def __init__(
        self,
        adata: AnnData | None = None,
        registry: dict | None = None,
        n_latent: int = 32,
        encoder_dims: Sequence[int] = (128, 128),
        decoder_dims: Sequence[int] = (128, 128),
        prior: Literal["normal"] = "normal",
        categorical_embedding_dims: dict[str, int] | None = None,
        **model_kwargs,
    ) -> None:
        super().__init__(adata, registry)

        for key in ["categorical_covariates", "batch_key"]:
            if key in model_kwargs:
                raise RuntimeError(
                    f"Passing {key} to DRVI model is no longer possible."
                    "It is enough to pass this argument to DRVI.setup_anndata."
                )

        n_batch = self.summary_stats.n_batch
        n_cats_per_cov = [n_batch] + self._get_n_cats_per_cov()
        n_continuous_cov = self.summary_stats.get("n_extra_continuous_covs", 0)

        categorical_covariates_dims = self._compute_categorical_covariates_dims(
            n_cats_per_cov, categorical_embedding_dims
        )

        self._module_kwargs = dict(
            n_input=self.summary_stats["n_vars"],
            n_latent=n_latent,
            encoder_dims=encoder_dims,
            decoder_dims=decoder_dims,
            n_cats_per_cov=n_cats_per_cov,
            n_continuous_cov=n_continuous_cov,
            prior=prior,
            categorical_covariate_dims=categorical_covariates_dims,
            **model_kwargs,
        )

        self.module = self._module_cls(**self._module_kwargs)

        self._model_summary_string = (
            "DRVI \n"
            f"Latent size: {self.module.n_latent}, "
            f"splits: {self.module.n_split_latent}, "
            f"pooling of splits: '{self.module.split_aggregation}', \n"
            f"Encoder dims: {encoder_dims}, \n"
            f"Decoder dims: {decoder_dims}, \n"
            f"Gene likelihood: {self.module.gene_likelihood}, \n"
        )

        # necessary line to get params that will be used for saving/loading
        self.init_params_ = self._get_init_params(locals())

        logger.info("The model has been initialized")

    def _get_n_cats_per_cov(self) -> list[int]:
        """Get the number of categories per categorical covariate."""
        if REGISTRY_KEYS.CAT_COVS_KEY in self.registry["field_registries"]:
            cat_cov_stats = self.registry["field_registries"][REGISTRY_KEYS.CAT_COVS_KEY]["state_registry"]
            return cat_cov_stats.get("n_cats_per_key", [])
        return []

    def _compute_categorical_covariates_dims(
        self,
        n_cats_per_cov: list[int],
        categorical_embedding_dims: dict[str, int] | None,
    ) -> list[int]:
        """Compute categorical covariates embedding dimensions.

        Parameters
        ----------
        n_cats_per_cov
            Number of categories per categorical covariate. This includes the batch covariate.
        categorical_embedding_dims
            Dictionary mapping covariate names to their embedding dimensions. This includes the batch covariate.

        Returns
        -------
        list[int]
            List of embedding dimensions for each categorical covariate.
        """
        categorical_embedding_dims = categorical_embedding_dims or {}
        categorical_covariates_dims = [self._DEFAULT_CATEGORICAL_EMBEDDING_DIM] * len(n_cats_per_cov)
        if n_cats_per_cov[0] > 1:
            batch_original_key = self.adata_manager.get_state_registry(REGISTRY_KEYS.BATCH_KEY).original_key
            if batch_original_key in categorical_embedding_dims:
                categorical_covariates_dims[0] = categorical_embedding_dims[batch_original_key]
        if len(n_cats_per_cov) > 1:
            cat_cov_names = self.adata_manager.get_state_registry(REGISTRY_KEYS.CAT_COVS_KEY).field_keys
            for i, obs_name in enumerate(cat_cov_names):
                if obs_name in categorical_embedding_dims:
                    categorical_covariates_dims[i + 1] = categorical_embedding_dims[obs_name]
        return categorical_covariates_dims

    @staticmethod
    def _update_source_registry_for_existing_model(source_registry: dict[str, Any]) -> dict[str, Any]:
        """Update the source registry for an existing model to the latest version if any updates are needed."""
        from packaging.version import Version

        source_registry_drvi_version = Version(
            source_registry.get("drvi_version", "0.1.0")
        )  # "0.1.0" for legacy code before pypi release
        logger.info(f"The model is trained with DRVI version {source_registry_drvi_version}.")

        while source_registry_drvi_version < Version(SCVI_VERSION):
            if source_registry_drvi_version < Version("0.1.10"):
                # No braking change up to 0.1.10
                source_registry_drvi_version = Version("0.1.10")
            elif source_registry_drvi_version == Version("0.1.10"):
                # log the transfer
                logger.info("Modifying model args from 0.1.10 to 0.1.11 (no user action required)")
                logger.info("Adding empty batch key ...")
                source_registry["setup_args"]["batch_key"] = None
                source_registry["field_registries"]["batch"] = {
                    "data_registry": {"attr_name": "obs", "attr_key": "_scvi_batch"},
                    "state_registry": {"categorical_mapping": np.array([0]), "original_key": "_scvi_batch"},
                    "summary_stats": {"n_batch": 1},
                }
                source_registry_drvi_version = Version("0.1.11")
                logger.info("Done updating source registry from 0.1.10 to 0.1.11.")
            else:
                # No braking change yet!
                source_registry_drvi_version = Version(SCVI_VERSION)
        logger.info(f"Loading source in DRVI version {SCVI_VERSION}.")

        return source_registry

    @classmethod
    @setup_anndata_dsp.dedent
    def setup_anndata(
        cls,
        adata: AnnData,
        layer: str | None = None,
        is_count_data: bool = True,
        batch_key: str | None = None,
        labels_key: str | None = None,
        categorical_covariate_keys: list[str] | None = None,
        continuous_covariate_keys: list[str] | None = None,
        **kwargs,
    ) -> None:
        """
        %(summary)s.

        Parameters
        ----------
        %(param_adata)s
        %(param_labels_key)s
        %(param_layer)s
        %(param_batch_key)s
        %(param_cat_cov_keys)s
        %(param_cont_cov_keys)s

        Returns
        -------
        %(returns)s
        """
        setup_method_args = cls._get_setup_method_args(**locals())
        setup_method_args["drvi_version"] = SCVI_VERSION

        # Manipulate kwargs in case of version updates (only when loading a model).
        if "source_registry" in kwargs:
            kwargs["source_registry"] = cls._update_source_registry_for_existing_model(kwargs["source_registry"])

        anndata_fields = [
            LayerField(REGISTRY_KEYS.X_KEY, layer, is_count_data=is_count_data),
            CategoricalObsField(REGISTRY_KEYS.BATCH_KEY, batch_key),
            CategoricalObsField(REGISTRY_KEYS.LABELS_KEY, labels_key),
            CategoricalJointObsField(REGISTRY_KEYS.CAT_COVS_KEY, categorical_covariate_keys),
            NumericalJointObsField(REGISTRY_KEYS.CONT_COVS_KEY, continuous_covariate_keys),
        ]
        adata_manager = AnnDataManager(fields=anndata_fields, setup_method_args=setup_method_args)
        adata_manager.register_fields(adata, **kwargs)
        cls.register_manager(adata_manager)
