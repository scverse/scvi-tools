from typing import Optional, Tuple

import numpy as np
import pandas as pd
from anndata import AnnData

from scvi._compat import Literal
from scvi.data import register_tensor_from_anndata
from scvi.dataloaders import AnnDataLoader
from scvi.external.stereoscope._module import RNADeconv, SpatialDeconv
from scvi.lightning import TrainingPlan
from scvi.model.base import BaseModelClass


class RNAStereoscope(BaseModelClass):
    """
    Reimplementation of Stereoscope [Andersson20]_ for deconvolution of spatial transcriptomics from single-cell transcriptomics.

    https://github.com/almaan/stereoscope.

    Parameters
    ----------
    sc_adata
        single-cell AnnData object that has been registered via :func:`~scvi.data.setup_anndata`.
    use_gpu
        Use the GPU or not.
    **model_kwargs
        Keyword args for :class:`~scvi.external.RNADeconv`

    Examples
    --------
    >>> sc_adata = anndata.read_h5ad(path_to_sc_anndata)
    >>> scvi.data.setup_anndata(sc_adata, label_key="labels")
    >>> stereo = scvi.external.RNAStereoscope(sc_adata)
    >>> stereo.train()
    """

    def __init__(
        self,
        sc_adata: AnnData,
        use_gpu: bool = True,
        **model_kwargs,
    ):
        super(RNAStereoscope, self).__init__(sc_adata, use_gpu=use_gpu)
        self.n_genes = self.summary_stats["n_vars"]
        self.n_labels = self.summary_stats["n_labels"]
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
        use_gpu: Optional[bool] = None,
        train_size: float = 1,
        validation_size: Optional[float] = None,
        batch_size: int = 128,
        plan_kwargs: Optional[dict] = None,
        **kwargs,
    ):
        """
        Trains the model using amortized variational inference.

        Parameters
        ----------
        max_epochs
            Number of epochs to train for
        lr
            Learning rate for optimization.
        use_gpu
            If `True`, use the GPU if available.
        train_size
            Size of training set in the range [0.0, 1.0].
        validation_size
            Size of the test set. If `None`, defaults to 1 - `train_size`. If
            `train_size + validation_size < 1`, the remaining cells belong to a test set.
        batch_size
            Minibatch size to use during training.
        plan_kwargs
            Keyword args for :class:`~scvi.lightning.TrainingPlan`. Keyword arguments passed to
            `train()` will overwrite values present in `plan_kwargs`, when appropriate.
        **kwargs
            Other keyword args for :class:`~scvi.lightning.Trainer`.
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

    @property
    def _plan_class(self):
        return TrainingPlan

    @property
    def _data_loader_cls(self):
        return AnnDataLoader


class SpatialStereoscope(BaseModelClass):
    """
    Reimplementation of Stereoscope [Andersson20]_ for deconvolution of spatial transcriptomics from single-cell transcriptomics.

    https://github.com/almaan/stereoscope.

    Parameters
    ----------
    st_adata
        spatial transcriptomics AnnData object that has been registered via :func:`~scvi.data.setup_anndata`.
    sc_params
        parameters of the model learned from the single-cell RNA seq data for deconvolution.
    cell_type_mapping
        numpy array mapping for the cell types used in the deconvolution
    use_gpu
        Use the GPU or not.
    prior_weight
        how to reweight the minibatches for stochastic optimization. "n_obs" is the valid
        procedure, "minibatch" is the procedure implemented in Stereoscope.
    **model_kwargs
        Keyword args for :class:`~scvi.external.SpatialDeconv`

    Examples
    --------
    >>> st_adata = anndata.read_h5ad(path_to_st_anndata)
    >>> scvi.data.setup_anndata(st_adata)
    >>> st_adata.obs["indices"] = np.arange(st_adata.n_obs)
    >>> register_tensor_from_anndata(st_adata, "ind_x", "obs", "indices")
    >>> stereo = scvi.external.SpatialStereoscope(st_adata, sc_params, cell_type_mapping)
    >>> stereo.train()
    >>> st_adata.obs["deconv"] = stereo.get_proportions()
    """

    def __init__(
        self,
        st_adata: AnnData,
        sc_params: Tuple[np.ndarray],
        cell_type_mapping: np.ndarray,
        use_gpu: bool = True,
        prior_weight: Literal["n_obs", "minibatch"] = "n_obs",
        **model_kwargs,
    ):
        st_adata.obs["_indices"] = np.arange(st_adata.n_obs)
        register_tensor_from_anndata(st_adata, "ind_x", "obs", "_indices")
        super().__init__(st_adata, use_gpu=use_gpu)
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
        use_gpu: bool = True,
        prior_weight: Literal["n_obs", "minibatch"] = "n_obs",
        **model_kwargs,
    ):
        """
        Alternate constructor for exploiting a pre-trained model on RNA-seq data

        Parameters
        -----------
        st_adata
            registed anndata object
        sc_model
            trained RNADeconv model
        use_gpu
            Use the GPU or not.
        prior_weight
            how to reweight the minibatches for stochastic optimization. "n_obs" is the valid
            procedure, "minibatch" is the procedure implemented in Stereoscope.
        **model_kwargs
            Keyword args for :class:`~scvi.external.SpatialDeconv`
        """
        return cls(
            st_adata,
            sc_model.module.get_params(),
            sc_model.scvi_setup_dict_["categorical_mappings"]["_scvi_labels"][
                "mapping"
            ],
            use_gpu=use_gpu,
            prior_weight=prior_weight,
            **model_kwargs,
        )

    def get_proportions(self, keep_noise=False) -> np.ndarray:
        """
        Returns the estimated cell type proportion for the spatial data. Shape is n_cells x n_labels OR n_cells x (n_labels + 1) if keep_noise

        Parameters:
        -----------
        keep_noise
            whether to account for the noise term as a standalone cell type in the proportion estimate.
        """
        column_names = self.cell_type_mapping
        if keep_noise:
            column_names = column_names.append("noise_term")
        return pd.DataFrame(
            data=self.module.get_proportions(keep_noise),
            columns=column_names,
            index=self.adata.obs.index,
        )

    def train(
        self,
        max_epochs: int = 400,
        lr: float = 0.01,
        use_gpu: Optional[bool] = None,
        batch_size: int = 128,
        plan_kwargs: Optional[dict] = None,
        **kwargs,
    ):
        """
        Trains the model using amortized variational inference.

        Parameters
        ----------
        max_epochs
            Number of epochs to train for
        lr
            Learning rate for optimization.
        use_gpu
            If `True`, use the GPU if available.
        batch_size
            Minibatch size to use during training.
        plan_kwargs
            Keyword args for :class:`~scvi.lightning.TrainingPlan`. Keyword arguments passed to
            `train()` will overwrite values present in `plan_kwargs`, when appropriate.
        **kwargs
            Other keyword args for :class:`~scvi.lightning.Trainer`.
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

    @property
    def _plan_class(self):
        return TrainingPlan

    @property
    def _data_loader_cls(self):
        return AnnDataLoader
