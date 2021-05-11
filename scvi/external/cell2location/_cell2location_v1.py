import numpy as np
import pandas as pd
from anndata import AnnData

import scvi
from scvi.model.base import BaseModelClass

from ._base import Cell2locationBaseModule, PltExportMixin, TrainSampleMixin
from ._cell2location_v1_module import LocationModelLinearDependentWMultiExperimentModel


def intersect_var(adata, cell_state_df):
    """
    Subset adata and cell_state_df to common variables (rows in cell_state_df).

    Parameters
    ----------
    adata
        anndata object
    cell_state_df
        pd.DataFrame with variables (genes) in rows and cell types/factors/covariates in columns.

    Returns
    -------

    """
    intersect = np.intersect1d(cell_state_df.index, adata.var_names)
    return adata[:, intersect].copy(), cell_state_df.loc[intersect, :].copy()


class Cell2location(TrainSampleMixin, BaseModelClass, PltExportMixin):
    """
    Reimplementation of cell2location [Kleshchevnikov20]_ model. User-end model class.

    https://github.com/BayraktarLab/cell2location

    Parameters
    ----------
    adata
        spatial AnnData object that has been registered via :func:`~scvi.data.setup_anndata`.
    cell_state_df
        pd.DataFrame with reference expression signatures
    use_gpu
        Use the GPU?
    **model_kwargs
        Keyword args for :class:`~scvi.external.LocationModelLinearDependentWMultiExperimentModel`

    Examples
    --------
    TODO add example
    >>>
    """

    def __init__(
        self,
        adata: AnnData,
        cell_state_df: pd.DataFrame,
        batch_size=None,
        model=None,
        **model_kwargs,
    ):
        # add index for each cell (provided to pyro plate for correct minibatching)
        adata.obs["_indices"] = np.arange(adata.n_obs).astype("int64")
        scvi.data.register_tensor_from_anndata(
            adata,
            registry_key="ind_x",
            adata_attr_name="obs",
            adata_key_name="_indices",
        )

        super().__init__(adata)

        if model is None:
            model = LocationModelLinearDependentWMultiExperimentModel

        self.cell_state_df_ = cell_state_df
        self.n_factors_ = cell_state_df.shape[1]
        self.factor_names_ = cell_state_df.columns.values

        self.module = Cell2locationBaseModule(
            model=model,
            n_obs=self.summary_stats["n_cells"],
            n_vars=self.summary_stats["n_vars"],
            n_factors=self.n_factors_,
            n_batch=self.summary_stats["n_batch"],
            batch_size=batch_size,
            cell_state_mat=self.cell_state_df_.values.astype("float32"),
            **model_kwargs,
        )
        self._model_summary_string = f'scVI-cell2location Model with the following params: \nn_factors: {self.n_factors_} \nn_batch: {self.summary_stats["n_batch"]} '
        self.init_params_ = self._get_init_params(locals())
