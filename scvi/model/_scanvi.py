import logging
import warnings
from copy import deepcopy
from typing import List, Optional, Sequence, Union

import numpy as np
import pandas as pd
import torch
from anndata import AnnData

from scvi import REGISTRY_KEYS
from scvi._compat import Literal
from scvi.data import AnnDataManager
from scvi.data._constants import _SETUP_ARGS_KEY
from scvi.data._utils import get_anndata_attribute
from scvi.data.fields import (
    CategoricalJointObsField,
    CategoricalObsField,
    LabelsWithUnlabeledObsField,
    LayerField,
    NumericalJointObsField,
    NumericalObsField,
)
from scvi.model._utils import _init_library_size
from scvi.model.base import SemiSupervisedTrainingMixin
from scvi.module import SCANVAE
from scvi.utils import setup_anndata_dsp

from ._scvi import SCVI
from .base import ArchesMixin, BaseModelClass, RNASeqMixin, VAEMixin

logger = logging.getLogger(__name__)


class SCANVI(
    RNASeqMixin, VAEMixin, ArchesMixin, SemiSupervisedTrainingMixin, BaseModelClass
):
    """
    Single-cell annotation using variational inference [Xu21]_.

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

    1. :doc:`/tutorials/notebooks/harmonization`
    2. :doc:`/tutorials/notebooks/scarches_scvi_tools`
    3. :doc:`/tutorials/notebooks/seed_labeling`
    """

    def __init__(
        self,
        adata: AnnData,
        n_hidden: int = 128,
        n_latent: int = 10,
        n_layers: int = 1,
        dropout_rate: float = 0.1,
        dispersion: Literal["gene", "gene-batch", "gene-label", "gene-cell"] = "gene",
        gene_likelihood: Literal["zinb", "nb", "poisson"] = "zinb",
        **model_kwargs,
    ):
        super().__init__(adata)
        scanvae_model_kwargs = dict(model_kwargs)

        self._set_indices_and_labels()

        # ignores unlabeled catgegory
        n_labels = self.summary_stats.n_labels - 1
        n_cats_per_cov = (
            self.adata_manager.get_state_registry(
                REGISTRY_KEYS.CAT_COVS_KEY
            ).n_cats_per_key
            if REGISTRY_KEYS.CAT_COVS_KEY in self.adata_manager.data_registry
            else None
        )

        n_batch = self.summary_stats.n_batch
        use_size_factor_key = (
            REGISTRY_KEYS.SIZE_FACTOR_KEY in self.adata_manager.data_registry
        )
        library_log_means, library_log_vars = None, None
        if not use_size_factor_key:
            library_log_means, library_log_vars = _init_library_size(
                self.adata_manager, n_batch
            )

        self.module = SCANVAE(
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
            library_log_means=library_log_means,
            library_log_vars=library_log_vars,
            **scanvae_model_kwargs,
        )

        self.unsupervised_history_ = None
        self.semisupervised_history_ = None

        self._model_summary_string = (
            "ScanVI Model with the following params: \nunlabeled_category: {}, n_hidden: {}, n_latent: {}"
            ", n_layers: {}, dropout_rate: {}, dispersion: {}, gene_likelihood: {}"
        ).format(
            self.unlabeled_category_,
            n_hidden,
            n_latent,
            n_layers,
            dropout_rate,
            dispersion,
            gene_likelihood,
        )
        self.init_params_ = self._get_init_params(locals())
        self.was_pretrained = False

    @classmethod
    def from_scvi_model(
        cls,
        scvi_model: SCVI,
        unlabeled_category: str,
        labels_key: Optional[str] = None,
        adata: Optional[AnnData] = None,
        **scanvi_kwargs,
    ):
        """
        Initialize scanVI model with weights from pretrained :class:`~scvi.model.SCVI` model.

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
        scanvi_kwargs
            kwargs for scANVI model
        """
        scvi_model._check_if_trained(
            message="Passed in scvi model hasn't been trained yet."
        )

        scanvi_kwargs = dict(scanvi_kwargs)
        init_params = scvi_model.init_params_
        non_kwargs = init_params["non_kwargs"]
        kwargs = init_params["kwargs"]
        kwargs = {k: v for (i, j) in kwargs.items() for (k, v) in j.items()}
        for k, v in {**non_kwargs, **kwargs}.items():
            if k in scanvi_kwargs.keys():
                warnings.warn(
                    "Ignoring param '{}' as it was already passed in to ".format(k)
                    + "pretrained scvi model with value {}.".format(v)
                )
                del scanvi_kwargs[k]

        if adata is None:
            adata = scvi_model.adata
        else:
            # validate new anndata against old model
            scvi_model._validate_anndata(adata)

        scvi_setup_args = deepcopy(scvi_model.adata_manager.registry[_SETUP_ARGS_KEY])
        scvi_labels_key = scvi_setup_args["labels_key"]
        if labels_key is None and scvi_labels_key is None:
            raise ValueError(
                "A `labels_key` is necessary as the SCVI model was initialized without one."
            )
        if scvi_labels_key is None:
            scvi_setup_args.update(dict(labels_key=labels_key))
        cls.setup_anndata(
            adata,
            unlabeled_category=unlabeled_category,
            **scvi_setup_args,
        )
        scanvi_model = cls(adata, **non_kwargs, **kwargs, **scanvi_kwargs)
        scvi_state_dict = scvi_model.module.state_dict()
        scanvi_model.module.load_state_dict(scvi_state_dict, strict=False)
        scanvi_model.was_pretrained = True

        return scanvi_model

    def _set_indices_and_labels(self):
        """
        Set indices for labeled and unlabeled cells.
        """
        labels_state_registry = self.adata_manager.get_state_registry(
            REGISTRY_KEYS.LABELS_KEY
        )
        self.original_label_key = labels_state_registry.original_key
        self.unlabeled_category_ = labels_state_registry.unlabeled_category

        labels = get_anndata_attribute(
            self.adata,
            self.adata_manager.data_registry.labels.attr_name,
            self.original_label_key,
        ).ravel()
        self._label_mapping = labels_state_registry.categorical_mapping

        # set unlabeled and labeled indices
        self._unlabeled_indices = np.argwhere(
            labels == self.unlabeled_category_
        ).ravel()
        self._labeled_indices = np.argwhere(labels != self.unlabeled_category_).ravel()
        self._code_to_label = {i: l for i, l in enumerate(self._label_mapping)}

    def predict(
        self,
        adata: Optional[AnnData] = None,
        indices: Optional[Sequence[int]] = None,
        soft: bool = False,
        batch_size: Optional[int] = None,
    ) -> Union[np.ndarray, pd.DataFrame]:
        """
        Return cell label predictions.

        Parameters
        ----------
        adata
            AnnData object that has been registered via :meth:`~scvi.model.SCANVI.setup_anndata`.
        indices
            Return probabilities for each class label.
        soft
            If True, returns per class probabilities
        batch_size
            Minibatch size for data loading into model. Defaults to `scvi.settings.batch_size`.
        """
        adata = self._validate_anndata(adata)

        if indices is None:
            indices = np.arange(adata.n_obs)

        scdl = self._make_data_loader(
            adata=adata,
            indices=indices,
            batch_size=batch_size,
        )
        y_pred = []
        for _, tensors in enumerate(scdl):
            x = tensors[REGISTRY_KEYS.X_KEY]
            batch = tensors[REGISTRY_KEYS.BATCH_KEY]

            cont_key = REGISTRY_KEYS.CONT_COVS_KEY
            cont_covs = tensors[cont_key] if cont_key in tensors.keys() else None

            cat_key = REGISTRY_KEYS.CAT_COVS_KEY
            cat_covs = tensors[cat_key] if cat_key in tensors.keys() else None

            pred = self.module.classify(
                x, batch_index=batch, cat_covs=cat_covs, cont_covs=cont_covs
            )
            if not soft:
                pred = pred.argmax(dim=1)
            y_pred.append(pred.detach().cpu())

        y_pred = torch.cat(y_pred).numpy()
        if not soft:
            predictions = []
            for p in y_pred:
                predictions.append(self._code_to_label[p])

            return np.array(predictions)
        else:
            n_labels = len(pred[0])
            pred = pd.DataFrame(
                y_pred,
                columns=self._label_mapping[:n_labels],
                index=adata.obs_names[indices],
            )
            return pred

    @classmethod
    @setup_anndata_dsp.dedent
    def setup_anndata(
        cls,
        adata: AnnData,
        labels_key: str,
        unlabeled_category: Union[str, int, float],
        layer: Optional[str] = None,
        batch_key: Optional[str] = None,
        size_factor_key: Optional[str] = None,
        categorical_covariate_keys: Optional[List[str]] = None,
        continuous_covariate_keys: Optional[List[str]] = None,
        **kwargs,
    ):
        """
        %(summary)s.

        Parameters
        ----------
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
            LabelsWithUnlabeledObsField(
                REGISTRY_KEYS.LABELS_KEY, labels_key, unlabeled_category
            ),
            NumericalObsField(
                REGISTRY_KEYS.SIZE_FACTOR_KEY, size_factor_key, required=False
            ),
            CategoricalJointObsField(
                REGISTRY_KEYS.CAT_COVS_KEY, categorical_covariate_keys
            ),
            NumericalJointObsField(
                REGISTRY_KEYS.CONT_COVS_KEY, continuous_covariate_keys
            ),
        ]
        adata_manager = AnnDataManager(
            fields=anndata_fields, setup_method_args=setup_method_args
        )
        adata_manager.register_fields(adata, **kwargs)
        cls.register_manager(adata_manager)
