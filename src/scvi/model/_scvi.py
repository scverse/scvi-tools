from __future__ import annotations

import logging
import warnings
from typing import TYPE_CHECKING

import numpy as np

import scvi
from scvi import REGISTRY_KEYS, settings
from scvi.data import AnnDataManager
from scvi.data._constants import ADATA_MINIFY_TYPE
from scvi.data._utils import _get_adata_minify_type
from scvi.data.fields import (
    CategoricalJointObsField,
    CategoricalObsField,
    LayerField,
    NumericalJointObsField,
    NumericalObsField,
)
from scvi.model._utils import _init_library_size
from scvi.model.base import EmbeddingMixin, UnsupervisedTrainingMixin
from scvi.module import VAE
from scvi.utils import setup_anndata_dsp

from .base import ArchesMixin, BaseMinifiedModeModelClass, RNASeqMixin, VAEMixin

if TYPE_CHECKING:
    from typing import Literal

    from anndata import AnnData
    from lightning import LightningDataModule


_SCVI_LATENT_QZM = "_scvi_latent_qzm"
_SCVI_LATENT_QZV = "_scvi_latent_qzv"
_SCVI_OBSERVED_LIB_SIZE = "_scvi_observed_lib_size"

logger = logging.getLogger(__name__)


class SCVI(
    EmbeddingMixin,
    RNASeqMixin,
    VAEMixin,
    ArchesMixin,
    UnsupervisedTrainingMixin,
    BaseMinifiedModeModelClass,
):
    """single-cell Variational Inference :cite:p:`Lopez18`.

    Parameters
    ----------
    adata
        AnnData object that has been registered via :meth:`~scvi.model.SCVI.setup_anndata`. If
        ``None``, then the underlying module will not be initialized until training, and a
        :class:`~lightning.pytorch.core.LightningDataModule` must be passed in during training
        (``EXPERIMENTAL``).
    n_hidden
        Number of nodes per hidden layer.
    n_latent
        Dimensionality of the latent space.
    n_layers
        Number of hidden layers used for encoder and decoder NNs.
    dropout_rate
        Dropout rate for neural networks.
    dispersion
        One of the following:

        * ``'gene'`` - dispersion parameter of NB is constant per gene across cells
        * ``'gene-batch'`` - dispersion can differ between different batches
        * ``'gene-label'`` - dispersion can differ between different labels
        * ``'gene-cell'`` - dispersion can differ for every gene in every cell
    gene_likelihood
        One of:

        * ``'nb'`` - Negative binomial distribution
        * ``'zinb'`` - Zero-inflated negative binomial distribution
        * ``'poisson'`` - Poisson distribution
        * ``'normal'`` - ``EXPERIMENTAL`` Normal distribution
    latent_distribution
        One of:

        * ``'normal'`` - Normal distribution
        * ``'ln'`` - Logistic normal distribution (Normal(0, I) transformed by softmax)
    **kwargs
        Additional keyword arguments for :class:`~scvi.module.VAE`.

    Examples
    --------
    >>> adata = anndata.read_h5ad(path_to_anndata)
    >>> scvi.model.SCVI.setup_anndata(adata, batch_key="batch")
    >>> vae = scvi.model.SCVI(adata)
    >>> vae.train()
    >>> adata.obsm["X_scVI"] = vae.get_latent_representation()
    >>> adata.obsm["X_normalized_scVI"] = vae.get_normalized_expression()

    Notes
    -----
    See further usage examples in the following tutorials:

    1. :doc:`/tutorials/notebooks/quick_start/api_overview`
    2. :doc:`/tutorials/notebooks/scrna/harmonization`
    3. :doc:`/tutorials/notebooks/scrna/scarches_scvi_tools`
    4. :doc:`/tutorials/notebooks/scrna/scvi_in_R`

    See Also
    --------
    :class:`~scvi.module.VAE`
    """

    _module_cls = VAE
    _LATENT_QZM_KEY = "scvi_latent_qzm"
    _LATENT_QZV_KEY = "scvi_latent_qzv"

    def __init__(
        self,
        adata: AnnData | None = None,
        registry: dict | None = None,
        n_hidden: int = 128,
        n_latent: int = 10,
        n_layers: int = 1,
        dropout_rate: float = 0.1,
        dispersion: Literal["gene", "gene-batch", "gene-label", "gene-cell"] = "gene",
        gene_likelihood: Literal["zinb", "nb", "poisson", "normal"] = "zinb",
        latent_distribution: Literal["normal", "ln"] = "normal",
        datamodule: LightningDataModule | None = None,
        **kwargs,
    ):
        super().__init__(adata, registry)

        self._module_kwargs = {
            "n_hidden": n_hidden,
            "n_latent": n_latent,
            "n_layers": n_layers,
            "dropout_rate": dropout_rate,
            "dispersion": dispersion,
            "gene_likelihood": gene_likelihood,
            "latent_distribution": latent_distribution,
            **kwargs,
        }
        self._model_summary_string = (
            "SCVI model with the following parameters: \n"
            f"n_hidden: {n_hidden}, n_latent: {n_latent}, n_layers: {n_layers}, "
            f"dropout_rate: {dropout_rate}, dispersion: {dispersion}, "
            f"gene_likelihood: {gene_likelihood}, latent_distribution: {latent_distribution}."
        )
        self._module_init_on_train = False

        if self._module_init_on_train:
            self.module = None
            warnings.warn(
                "Model was initialized without `adata`. The module will be initialized when "
                "calling `train`. This behavior is experimental and may change in the future.",
                UserWarning,
                stacklevel=settings.warnings_stacklevel,
            )
        else:
            if adata is not None:
                n_cats_per_cov = (
                    self.adata_manager.get_state_registry(
                        REGISTRY_KEYS.CAT_COVS_KEY
                    ).n_cats_per_key
                    if REGISTRY_KEYS.CAT_COVS_KEY in self.adata_manager.data_registry
                    else None
                )
            else:
                # custom datamodule
                n_cats_per_cov = self.summary_stats[f"n_{REGISTRY_KEYS.CAT_COVS_KEY}"]
                if n_cats_per_cov == 0:
                    n_cats_per_cov = None

            n_batch = self.summary_stats.n_batch
            use_size_factor_key = (
                False  # REGISTRY_KEYS.SIZE_FACTOR_KEY in self.adata_manager.data_registry
            )
            library_log_means, library_log_vars = None, None
            if (
                use_size_factor_key
                and self.minified_data_type != ADATA_MINIFY_TYPE.LATENT_POSTERIOR
            ):
                library_log_means, library_log_vars = _init_library_size(
                    self.adata_manager, n_batch
                )
            self.module = self._module_cls(
                n_input=self.summary_stats.n_vars,
                n_batch=n_batch,
                n_labels=self.summary_stats.n_labels,
                n_continuous_cov=self.summary_stats.get("n_extra_continuous_covs", 0),
                n_cats_per_cov=n_cats_per_cov,
                n_hidden=n_hidden,
                n_latent=n_latent,
                n_layers=n_layers,
                dropout_rate=dropout_rate,
                dispersion=dispersion,
                gene_likelihood=gene_likelihood,
                latent_distribution=latent_distribution,
                use_size_factor_key=use_size_factor_key,
                library_log_means=library_log_means,
                library_log_vars=library_log_vars,
                **kwargs,
            )
            self.module.minified_data_type = self.minified_data_type

        self.init_params_ = self._get_init_params(locals())

    @classmethod
    @setup_anndata_dsp.dedent
    def setup_anndata(
        cls,
        adata: AnnData,
        layer: str | None = None,
        batch_key: str | None = None,
        labels_key: str | None = None,
        size_factor_key: str | None = None,
        categorical_covariate_keys: list[str] | None = None,
        continuous_covariate_keys: list[str] | None = None,
        **kwargs,
    ):
        """%(summary)s.

        Parameters
        ----------
        %(param_adata)s
        %(param_layer)s
        %(param_batch_key)s
        %(param_labels_key)s
        %(param_size_factor_key)s
        %(param_cat_cov_keys)s
        %(param_cont_cov_keys)s
        """
        setup_method_args = cls._get_setup_method_args(**locals())
        anndata_fields = [
            LayerField(REGISTRY_KEYS.X_KEY, layer, is_count_data=True),
            CategoricalObsField(REGISTRY_KEYS.BATCH_KEY, batch_key),
            CategoricalObsField(REGISTRY_KEYS.LABELS_KEY, labels_key),
            NumericalObsField(REGISTRY_KEYS.SIZE_FACTOR_KEY, size_factor_key, required=False),
            CategoricalJointObsField(REGISTRY_KEYS.CAT_COVS_KEY, categorical_covariate_keys),
            NumericalJointObsField(REGISTRY_KEYS.CONT_COVS_KEY, continuous_covariate_keys),
        ]
        # register new fields if the adata is minified
        adata_minify_type = _get_adata_minify_type(adata)
        if adata_minify_type is not None:
            anndata_fields += cls._get_fields_for_adata_minification(adata_minify_type)
        adata_manager = AnnDataManager(fields=anndata_fields, setup_method_args=setup_method_args)
        adata_manager.register_fields(adata, **kwargs)
        cls.register_manager(adata_manager)

    @classmethod
    @setup_anndata_dsp.dedent
    def setup_datamodule(
        cls,
        datamodule: LightningDataModule | None = None,
        source_registry=None,
        layer: str | None = None,
        batch_key: list[str] | None = None,
        labels_key: str | None = None,
        size_factor_key: str | None = None,
        categorical_covariate_keys: list[str] | None = None,
        continuous_covariate_keys: list[str] | None = None,
        **kwargs,
    ):
        """%(summary)s.

        Parameters
        ----------
        %(param_datamodule)s
        %(param_source_registry)s
        %(param_layer)s
        %(param_batch_key)s
        %(param_labels_key)s
        %(param_size_factor_key)s
        %(param_cat_cov_keys)s
        %(param_cont_cov_keys)s
        """
        if datamodule.__class__.__name__ == "CensusSCVIDataModule":
            # CZI
            categorical_mapping = datamodule.datapipe.obs_encoders["batch"].classes_
            column_names = list(
                datamodule.datapipe.var_query.coords[0]
                if datamodule.datapipe.var_query is not None
                else range(datamodule.n_vars)
            )
            n_batch = datamodule.n_batch
        else:
            # Anndata -> CZI
            # if we are here and datamodule is actually an AnnData object
            # it means we init the custom dataloder model with anndata
            categorical_mapping = source_registry["field_registries"]["batch"]["state_registry"][
                "categorical_mapping"
            ]
            column_names = datamodule.var.soma_joinid.values
            n_batch = source_registry["field_registries"]["batch"]["summary_stats"]["n_batch"]

        datamodule.registry = {
            "scvi_version": scvi.__version__,
            "model_name": "SCVI",
            "setup_args": {
                "layer": layer,
                "batch_key": batch_key,
                "labels_key": labels_key,
                "size_factor_key": size_factor_key,
                "categorical_covariate_keys": categorical_covariate_keys,
                "continuous_covariate_keys": continuous_covariate_keys,
            },
            "field_registries": {
                "X": {
                    "data_registry": {"attr_name": "X", "attr_key": None},
                    "state_registry": {
                        "n_obs": datamodule.n_obs,
                        "n_vars": datamodule.n_vars,
                        "column_names": [str(i) for i in column_names],
                    },
                    "summary_stats": {"n_vars": datamodule.n_vars, "n_cells": datamodule.n_obs},
                },
                "batch": {
                    "data_registry": {"attr_name": "obs", "attr_key": "_scvi_batch"},
                    "state_registry": {
                        "categorical_mapping": categorical_mapping,
                        "original_key": "batch",
                    },
                    "summary_stats": {"n_batch": n_batch},
                },
                "labels": {
                    "data_registry": {"attr_name": "obs", "attr_key": "_scvi_labels"},
                    "state_registry": {
                        "categorical_mapping": np.array([0]),
                        "original_key": "_scvi_labels",
                    },
                    "summary_stats": {"n_labels": 1},
                },
                "size_factor": {"data_registry": {}, "state_registry": {}, "summary_stats": {}},
                "extra_categorical_covs": {
                    "data_registry": {},
                    "state_registry": {},
                    "summary_stats": {"n_extra_categorical_covs": 0},
                },
                "extra_continuous_covs": {
                    "data_registry": {},
                    "state_registry": {},
                    "summary_stats": {"n_extra_continuous_covs": 0},
                },
            },
            "setup_method_name": "setup_datamodule",
        }
