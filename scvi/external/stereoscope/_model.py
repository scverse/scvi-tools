import logging
from typing import Literal, Optional, Tuple, Union

import numpy as np
import pandas as pd
import torch
from anndata import AnnData

from scvi import REGISTRY_KEYS
from scvi.data import AnnDataManager
from scvi.data.fields import CategoricalObsField, LayerField, NumericalObsField
from scvi.external.stereoscope._module import RNADeconv, SpatialDeconv
from scvi.model.base import BaseModelClass, UnsupervisedTrainingMixin
from scvi.utils import setup_anndata_dsp

logger = logging.getLogger(__name__)


class RNAStereoscope(UnsupervisedTrainingMixin, BaseModelClass):
    """
    Reimplementation of Stereoscope :cite:p:`Andersson20` for deconvolution of spatial transcriptomics from single-cell transcriptomics.

    https://github.com/almaan/stereoscope.

    Parameters
    ----------
    sc_adata
        single-cell AnnData object that has been registered via :meth:`~scvi.external.RNAStereoscope.setup_anndata`.
    **model_kwargs
        Keyword args for :class:`~scvi.external.stereoscope.RNADeconv`

    Examples
    --------
    >>> sc_adata = anndata.read_h5ad(path_to_sc_anndata)
    >>> scvi.external.RNAStereoscope.setup_anndata(sc_adata, labels_key="labels")
    >>> stereo = scvi.external.stereoscope.RNAStereoscope(sc_adata)
    >>> stereo.train()
    """

    def __init__(
        self,
        sc_adata: AnnData,
        **model_kwargs,
    ):
        super().__init__(sc_adata)
        self.n_genes = self.summary_stats.n_vars
        self.n_labels = self.summary_stats.n_labels
        # first we have the scRNA-seq model
        self.module = RNADeconv(
            n_genes=self.n_genes,
            n_labels=self.n_labels,
            **model_kwargs,
        )
        self._model_summary_string = (
            "RNADeconv Model with params: \nn_genes: {}, n_labels: {}"
        ).format(
            self.n_genes,
            self.n_labels,
        )
        self.init_params_ = self._get_init_params(locals())

    def train(
        self,
        max_epochs: int = 400,
        lr: float = 0.01,
        use_gpu: Optional[Union[str, int, bool]] = None,
        train_size: float = 1,
        validation_size: Optional[float] = None,
        batch_size: int = 128,
        plan_kwargs: Optional[dict] = None,
        **kwargs,
    ):
        """
        Trains the model using MAP inference.

        Parameters
        ----------
        max_epochs
            Number of epochs to train for
        lr
            Learning rate for optimization.
        use_gpu
            Use default GPU if available (if None or True), or index of GPU to use (if int),
            or name of GPU (if str, e.g., `'cuda:0'`), or use CPU (if False).
        train_size
            Size of training set in the range [0.0, 1.0].
        validation_size
            Size of the test set. If `None`, defaults to 1 - `train_size`. If
            `train_size + validation_size < 1`, the remaining cells belong to a test set.
        batch_size
            Minibatch size to use during training.
        plan_kwargs
            Keyword args for :class:`~scvi.train.TrainingPlan`. Keyword arguments passed to
            `train()` will overwrite values present in `plan_kwargs`, when appropriate.
        **kwargs
            Other keyword args for :class:`~scvi.train.Trainer`.
        """
        update_dict = {
            "lr": lr,
        }
        if plan_kwargs is not None:
            plan_kwargs.update(update_dict)
        else:
            plan_kwargs = update_dict
        super().train(
            max_epochs=max_epochs,
            use_gpu=use_gpu,
            train_size=train_size,
            validation_size=validation_size,
            batch_size=batch_size,
            plan_kwargs=plan_kwargs,
            **kwargs,
        )

    @classmethod
    @setup_anndata_dsp.dedent
    def setup_anndata(
        cls,
        adata: AnnData,
        labels_key: Optional[str] = None,
        layer: Optional[str] = None,
        **kwargs,
    ):
        """
        %(summary)s.

        Parameters
        ----------
        %(param_labels_key)s
        %(param_layer)s
        """
        setup_method_args = cls._get_setup_method_args(**locals())
        anndata_fields = [
            LayerField(REGISTRY_KEYS.X_KEY, layer, is_count_data=True),
            CategoricalObsField(REGISTRY_KEYS.LABELS_KEY, labels_key),
        ]
        adata_manager = AnnDataManager(
            fields=anndata_fields, setup_method_args=setup_method_args
        )
        adata_manager.register_fields(adata, **kwargs)
        cls.register_manager(adata_manager)


class SpatialStereoscope(UnsupervisedTrainingMixin, BaseModelClass):
    """
    Reimplementation of Stereoscope :cite:p:`Andersson20` for deconvolution of spatial transcriptomics from single-cell transcriptomics.

    https://github.com/almaan/stereoscope.

    Parameters
    ----------
    st_adata
        spatial transcriptomics AnnData object that has been registered via :meth:`~scvi.external.SpatialStereoscope.setup_anndata`.
    sc_params
        parameters of the model learned from the single-cell RNA seq data for deconvolution.
    cell_type_mapping
        numpy array mapping for the cell types used in the deconvolution
    prior_weight
        how to reweight the minibatches for stochastic optimization. "n_obs" is the valid
        procedure, "minibatch" is the procedure implemented in Stereoscope.
    **model_kwargs
        Keyword args for :class:`~scvi.external.stereoscope.SpatialDeconv`

    Examples
    --------
    >>> sc_adata = anndata.read_h5ad(path_to_sc_anndata)
    >>> scvi.external.RNAStereoscope.setup_anndata(sc_adata, labels_key="labels")
    >>> sc_model = scvi.external.stereoscope.RNAStereoscope(sc_adata)
    >>> sc_model.train()
    >>> st_adata = anndata.read_h5ad(path_to_st_anndata)
    >>> scvi.external.SpatialStereoscope.setup_anndata(st_adata)
    >>> stereo = scvi.external.stereoscope.SpatialStereoscope.from_rna_model(st_adata, sc_model)
    >>> stereo.train()
    >>> st_adata.obsm["deconv"] = stereo.get_proportions()

    Notes
    -----
    See further usage examples in the following tutorials:

    1. :doc:`/user_guide/notebooks/stereoscope_heart_LV_tutorial`
    """

    def __init__(
        self,
        st_adata: AnnData,
        sc_params: Tuple[np.ndarray],
        cell_type_mapping: np.ndarray,
        prior_weight: Literal["n_obs", "minibatch"] = "n_obs",
        **model_kwargs,
    ):
        super().__init__(st_adata)
        self.module = SpatialDeconv(
            n_spots=st_adata.n_obs,
            sc_params=sc_params,
            prior_weight=prior_weight,
            **model_kwargs,
        )
        self._model_summary_string = (
            "RNADeconv Model with params: \nn_spots: {}"
        ).format(
            st_adata.n_obs,
        )
        self.cell_type_mapping = cell_type_mapping
        self.init_params_ = self._get_init_params(locals())

    @classmethod
    def from_rna_model(
        cls,
        st_adata: AnnData,
        sc_model: RNAStereoscope,
        prior_weight: Literal["n_obs", "minibatch"] = "n_obs",
        layer: Optional[str] = None,
        **model_kwargs,
    ):
        """
        Alternate constructor for exploiting a pre-trained model on RNA-seq data.

        Parameters
        ----------
        st_adata
            registed anndata object
        sc_model
            trained RNADeconv model
        prior_weight
            how to reweight the minibatches for stochastic optimization. "n_obs" is the valid
            procedure, "minibatch" is the procedure implemented in Stereoscope.
        layer
            if not `None`, uses this as the key in `adata.layers` for raw count data.
        **model_kwargs
            Keyword args for :class:`~scvi.external.SpatialDeconv`
        """
        cls.setup_anndata(st_adata, layer=layer)
        return cls(
            st_adata,
            sc_model.module.get_params(),
            sc_model.adata_manager.get_state_registry(
                REGISTRY_KEYS.LABELS_KEY
            ).categorical_mapping,
            prior_weight=prior_weight,
            **model_kwargs,
        )

    def get_proportions(self, keep_noise=False) -> pd.DataFrame:
        """
        Returns the estimated cell type proportion for the spatial data.

        Shape is n_cells x n_labels OR n_cells x (n_labels + 1) if keep_noise

        Parameters
        ----------
        keep_noise
            whether to account for the noise term as a standalone cell type in the proportion estimate.
        """
        self._check_if_trained()

        column_names = self.cell_type_mapping
        if keep_noise:
            column_names = column_names.append("noise_term")
        return pd.DataFrame(
            data=self.module.get_proportions(keep_noise),
            columns=column_names,
            index=self.adata.obs.index,
        )

    def get_scale_for_ct(
        self,
        y: np.ndarray,
    ) -> np.ndarray:
        r"""
        Calculate the cell type specific expression.

        Parameters
        ----------
        y
            numpy array containing the list of cell types
        Returns
        -------
        gene_expression
        """
        self._check_if_trained()
        ind_y = np.array([np.where(ct == self.cell_type_mapping)[0][0] for ct in y])
        if ind_y.shape != y.shape:
            raise ValueError(
                "Incorrect shape after matching cell types to reference mapping. Please check cell type query."
            )
        px_scale = self.module.get_ct_specific_expression(torch.tensor(ind_y)[:, None])
        return np.array(px_scale.cpu())

    def train(
        self,
        max_epochs: int = 400,
        lr: float = 0.01,
        use_gpu: Optional[Union[str, int, bool]] = None,
        batch_size: int = 128,
        plan_kwargs: Optional[dict] = None,
        **kwargs,
    ):
        """
        Trains the model using MAP inference.

        Parameters
        ----------
        max_epochs
            Number of epochs to train for
        lr
            Learning rate for optimization.
        use_gpu
            Use default GPU if available (if None or True), or index of GPU to use (if int),
            or name of GPU (if str, e.g., `'cuda:0'`), or use CPU (if False).
        batch_size
            Minibatch size to use during training.
        plan_kwargs
            Keyword args for :class:`~scvi.train.TrainingPlan`. Keyword arguments passed to
            `train()` will overwrite values present in `plan_kwargs`, when appropriate.
        **kwargs
            Other keyword args for :class:`~scvi.train.Trainer`.
        """
        update_dict = {
            "lr": lr,
        }
        if plan_kwargs is not None:
            plan_kwargs.update(update_dict)
        else:
            plan_kwargs = update_dict
        super().train(
            max_epochs=max_epochs,
            use_gpu=use_gpu,
            train_size=1,
            validation_size=None,
            batch_size=batch_size,
            plan_kwargs=plan_kwargs,
            **kwargs,
        )

    @classmethod
    @setup_anndata_dsp.dedent
    def setup_anndata(
        cls,
        adata: AnnData,
        layer: Optional[str] = None,
        **kwargs,
    ):
        """
        %(summary)s.

        Parameters
        ----------
        %(param_layer)s
        """
        setup_method_args = cls._get_setup_method_args(**locals())
        # add index for each cell (provided to pyro plate for correct minibatching)
        adata.obs["_indices"] = np.arange(adata.n_obs)
        anndata_fields = [
            LayerField(REGISTRY_KEYS.X_KEY, layer, is_count_data=True),
            NumericalObsField(REGISTRY_KEYS.INDICES_KEY, "_indices"),
        ]
        adata_manager = AnnDataManager(
            fields=anndata_fields, setup_method_args=setup_method_args
        )
        adata_manager.register_fields(adata, **kwargs)
        cls.register_manager(adata_manager)
