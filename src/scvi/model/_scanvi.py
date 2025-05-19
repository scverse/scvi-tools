from __future__ import annotations

import logging
import warnings
from copy import deepcopy
from typing import TYPE_CHECKING

from scvi import REGISTRY_KEYS, settings
from scvi.data import AnnDataManager
from scvi.data._constants import (
    _SETUP_ARGS_KEY,
    ADATA_MINIFY_TYPE,
)
from scvi.data._utils import _get_adata_minify_type, _is_minified
from scvi.data.fields import (
    CategoricalJointObsField,
    CategoricalObsField,
    LabelsWithUnlabeledObsField,
    LayerField,
    NumericalJointObsField,
    NumericalObsField,
)
from scvi.model._utils import _init_library_size
from scvi.module import SCANVAE
from scvi.train import SemiSupervisedTrainingPlan
from scvi.utils import setup_anndata_dsp

from .base import (
    ArchesMixin,
    BaseMinifiedModeModelClass,
    RNASeqMixin,
    SemisupervisedTrainingMixin,
    VAEMixin,
)

if TYPE_CHECKING:
    from typing import Literal

    from anndata import AnnData
    from lightning import LightningDataModule

    from ._scvi import SCVI

_SCANVI_LATENT_QZM = "_scanvi_latent_qzm"
_SCANVI_LATENT_QZV = "_scanvi_latent_qzv"
_SCANVI_OBSERVED_LIB_SIZE = "_scanvi_observed_lib_size"

logger = logging.getLogger(__name__)


class SCANVI(
    RNASeqMixin, SemisupervisedTrainingMixin, VAEMixin, ArchesMixin, BaseMinifiedModeModelClass
):
    """Single-cell annotation using variational inference :cite:p:`Xu21`.

    Inspired from M1 + M2 model, as described in (https://arxiv.org/pdf/1406.5298.pdf).

    Parameters
    ----------
    adata
        AnnData object that has been registered via :meth:`~scvi.model.SCANVI.setup_anndata`.
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
    use_observed_lib_size
        If ``True``, use the observed library size for RNA as the scaling factor in the mean of the
        conditional distribution.
    linear_classifier
        If ``True``, uses a single linear layer for classification instead of a
        multi-layer perceptron.
    **model_kwargs
        Keyword args for :class:`~scvi.module.SCANVAE`

    Examples
    --------
    >>> adata = anndata.read_h5ad(path_to_anndata)
    >>> scvi.model.SCANVI.setup_anndata(adata, batch_key="batch", labels_key="labels")
    >>> vae = scvi.model.SCANVI(adata, "Unknown")
    >>> vae.train()
    >>> adata.obsm["X_scVI"] = vae.get_latent_representation()
    >>> adata.obs["pred_label"] = vae.predict()

    Notes
    -----
    See further usage examples in the following tutorials:

    1. :doc:`/tutorials/notebooks/scrna/harmonization`
    2. :doc:`/tutorials/notebooks/scrna/scarches_scvi_tools`
    3. :doc:`/tutorials/notebooks/scrna/seed_labeling`
    """

    _module_cls = SCANVAE
    _training_plan_cls = SemiSupervisedTrainingPlan
    _LATENT_QZM_KEY = "scanvi_latent_qzm"
    _LATENT_QZV_KEY = "scanvi_latent_qzv"

    def __init__(
        self,
        adata: AnnData | None = None,
        registry: dict | None = None,
        n_hidden: int = 128,
        n_latent: int = 10,
        n_layers: int = 1,
        dropout_rate: float = 0.1,
        dispersion: Literal["gene", "gene-batch", "gene-label", "gene-cell"] = "gene",
        gene_likelihood: Literal["zinb", "nb", "poisson"] = "zinb",
        use_observed_lib_size: bool = True,
        linear_classifier: bool = False,
        datamodule: LightningDataModule | None = None,
        **model_kwargs,
    ):
        super().__init__(adata, registry)
        scanvae_model_kwargs = dict(model_kwargs)

        self._set_indices_and_labels(datamodule)

        # ignores unlabeled category if inside the labels
        if self.unlabeled_category_ is not None and self.unlabeled_category_ in self.labels_:
            n_labels = self.summary_stats.n_labels - 1
        else:
            if adata is not None and len(set(self.labels_)) == (self.summary_stats.n_labels - 1):
                n_labels = self.summary_stats.n_labels - 1
            else:
                n_labels = self.summary_stats.n_labels
        if adata is not None:
            n_cats_per_cov = (
                self.adata_manager.get_state_registry(REGISTRY_KEYS.CAT_COVS_KEY).n_cats_per_key
                if REGISTRY_KEYS.CAT_COVS_KEY in self.adata_manager.data_registry
                else None
            )
        else:
            # custom datamodule
            if (
                len(
                    self.registry["field_registries"][f"{REGISTRY_KEYS.CAT_COVS_KEY}"][
                        "state_registry"
                    ]
                )
                > 0
            ):
                n_cats_per_cov = tuple(
                    self.registry["field_registries"][f"{REGISTRY_KEYS.CAT_COVS_KEY}"][
                        "state_registry"
                    ]["n_cats_per_key"]
                )
            else:
                n_cats_per_cov = None

        n_batch = self.summary_stats.n_batch
        use_size_factor_key = self.registry_["setup_args"][f"{REGISTRY_KEYS.SIZE_FACTOR_KEY}_key"]
        library_log_means, library_log_vars = None, None
        if (
            not use_size_factor_key
            and self.minified_data_type != ADATA_MINIFY_TYPE.LATENT_POSTERIOR
            and not use_observed_lib_size
        ):
            library_log_means, library_log_vars = _init_library_size(self.adata_manager, n_batch)

        self.module = self._module_cls(
            n_input=self.summary_stats.n_vars,
            n_batch=n_batch,
            n_labels=n_labels,
            n_continuous_cov=self.summary_stats.get("n_extra_continuous_covs", 0),
            n_cats_per_cov=n_cats_per_cov,
            n_hidden=n_hidden,
            n_latent=n_latent,
            n_layers=n_layers,
            dropout_rate=dropout_rate,
            dispersion=dispersion,
            gene_likelihood=gene_likelihood,
            use_size_factor_key=use_size_factor_key,
            use_observed_lib_size=use_observed_lib_size,
            library_log_means=library_log_means,
            library_log_vars=library_log_vars,
            linear_classifier=linear_classifier,
            **scanvae_model_kwargs,
        )
        self.module.minified_data_type = self.minified_data_type

        self.unsupervised_history_ = None
        self.semisupervised_history_ = None

        self._model_summary_string = (
            f"ScanVI Model with the following params: \n"
            f"unlabeled_category: {self.unlabeled_category_}, n_hidden: {n_hidden}, "
            f"n_latent: {n_latent}, n_layers: {n_layers}, dropout_rate: {dropout_rate}, "
            f"dispersion: {dispersion}, gene_likelihood: {gene_likelihood}"
        )
        self.init_params_ = self._get_init_params(locals())
        self.was_pretrained = False
        self.n_labels = n_labels

    @classmethod
    def from_scvi_model(
        cls,
        scvi_model: SCVI,
        unlabeled_category: str,
        labels_key: str | None = None,
        adata: AnnData | None = None,
        registry: dict | None = None,
        **scanvi_kwargs,
    ):
        """Initialize scanVI model with weights from pretrained :class:`~scvi.model.SCVI` model.

        Parameters
        ----------
        scvi_model
            Pretrained scvi model
        labels_key
            key in `adata.obs` for label information. Label categories can not be different if
            labels_key was used to setup the SCVI model. If None, uses the `labels_key` used to
            setup the SCVI model. If that was None, and error is raised.
        unlabeled_category
            Value used for unlabeled cells in `labels_key` used to setup AnnData with scvi.
        adata
            AnnData object that has been registered via :meth:`~scvi.model.SCANVI.setup_anndata`.
        registry
            Registry of the datamodule used to train scANVI model.
        scanvi_kwargs
            kwargs for scANVI model
        """
        scvi_model._check_if_trained(message="Passed in scvi model hasn't been trained yet.")

        scanvi_kwargs = dict(scanvi_kwargs)
        init_params = scvi_model.init_params_
        non_kwargs = init_params["non_kwargs"]
        kwargs = init_params["kwargs"]
        kwargs = {k: v for (i, j) in kwargs.items() for (k, v) in j.items()}
        for k, v in {**non_kwargs, **kwargs}.items():
            if k in scanvi_kwargs.keys():
                warnings.warn(
                    f"Ignoring param '{k}' as it was already passed in to pretrained "
                    f"SCVI model with value {v}.",
                    UserWarning,
                    stacklevel=settings.warnings_stacklevel,
                )
                del scanvi_kwargs[k]

        if scvi_model.minified_data_type == ADATA_MINIFY_TYPE.LATENT_POSTERIOR:
            raise ValueError(
                "We cannot use the given scVI model to initialize scANVI because it has "
                "minified adata. Keep counts when minifying model using "
                "minified_data_type='latent_posterior_parameters_with_counts'."
            )

        if adata is None:
            adata = scvi_model.adata
        elif adata:
            if _is_minified(adata):
                raise ValueError("Please provide a non-minified `adata` to initialize scANVI.")
            # validate new anndata against old model
            scvi_model._validate_anndata(adata)
        else:
            adata = None

        scvi_setup_args = deepcopy(scvi_model.registry[_SETUP_ARGS_KEY])
        scvi_labels_key = scvi_setup_args["labels_key"]
        if labels_key is None and scvi_labels_key is None:
            raise ValueError(
                "A `labels_key` is necessary as the scVI model was initialized without one."
            )
        if scvi_labels_key is None:
            scvi_setup_args.update({"labels_key": labels_key})
        if adata is not None:
            cls.setup_anndata(
                adata,
                unlabeled_category=unlabeled_category,
                use_minified=False,
                **scvi_setup_args,
            )

        scanvi_model = cls(adata, scvi_model.registry, **non_kwargs, **kwargs, **scanvi_kwargs)
        scvi_state_dict = scvi_model.module.state_dict()
        scanvi_model.module.load_state_dict(scvi_state_dict, strict=False)
        scanvi_model.was_pretrained = True

        return scanvi_model

    @classmethod
    @setup_anndata_dsp.dedent
    def setup_anndata(
        cls,
        adata: AnnData,
        labels_key: str,
        unlabeled_category: str,
        layer: str | None = None,
        batch_key: str | None = None,
        size_factor_key: str | None = None,
        categorical_covariate_keys: list[str] | None = None,
        continuous_covariate_keys: list[str] | None = None,
        use_minified: bool = True,
        **kwargs,
    ):
        """%(summary)s.

        Parameters
        ----------
        %(param_adata)s
        %(param_labels_key)s
        %(param_unlabeled_category)s
        %(param_layer)s
        %(param_batch_key)s
        %(param_size_factor_key)s
        %(param_cat_cov_keys)s
        %(param_cont_cov_keys)s
        use_minified
            If True, will register the minified version of the adata if possible.
        """
        setup_method_args = cls._get_setup_method_args(**locals())
        anndata_fields = [
            LayerField(REGISTRY_KEYS.X_KEY, layer, is_count_data=True),
            CategoricalObsField(REGISTRY_KEYS.BATCH_KEY, batch_key),
            LabelsWithUnlabeledObsField(REGISTRY_KEYS.LABELS_KEY, labels_key, unlabeled_category),
            NumericalObsField(REGISTRY_KEYS.SIZE_FACTOR_KEY, size_factor_key, required=False),
            CategoricalJointObsField(REGISTRY_KEYS.CAT_COVS_KEY, categorical_covariate_keys),
            NumericalJointObsField(REGISTRY_KEYS.CONT_COVS_KEY, continuous_covariate_keys),
        ]
        # register new fields if the adata is minified
        if adata:
            adata_minify_type = _get_adata_minify_type(adata)
            if adata_minify_type is not None and use_minified:
                anndata_fields += cls._get_fields_for_adata_minification(adata_minify_type)
            adata_manager = AnnDataManager(
                fields=anndata_fields, setup_method_args=setup_method_args
            )
            adata_manager.register_fields(adata, **kwargs)
            cls.register_manager(adata_manager)
