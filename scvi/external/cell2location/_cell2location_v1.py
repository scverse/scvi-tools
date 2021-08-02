from typing import Optional

import numpy as np
import pandas as pd
from anndata import AnnData
from pyro import clear_param_store
from pyro.nn import PyroModule

import scvi
from scvi import _CONSTANTS
from scvi.data._anndata import get_from_registry
from scvi.model.base import BaseModelClass, PyroSampleMixin, PyroSviTrainMixin

from ._base import PltExportMixin, QuantileMixin
from ._cell2location_v1_module import Cell2locationBaseModule
from ._cell2location_v3_module import (
    LocationModelLinearDependentWMultiExperimentLocationBackgroundNormLevelGeneAlphaPyroModel,
)


class Cell2location(
    QuantileMixin, PyroSampleMixin, PyroSviTrainMixin, PltExportMixin, BaseModelClass
):
    """
    Reimplementation of cell2location [Kleshchevnikov20]_ model. User-end model class.

    https://github.com/BayraktarLab/cell2location

    Parameters
    ----------
    adata
        spatial AnnData object that has been registered via :func:`~scvi.data.setup_anndata`.
    cell_state_df
        pd.DataFrame with reference expression signatures for each gene (rows) in each cell type/population (columns).
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
        model_class: Optional[PyroModule] = None,
        detection_mean_per_sample: bool = False,
        detection_mean_correction: float = 1.0,
        **model_kwargs,
    ):
        # in case any other model was created before that shares the same parameter names.
        clear_param_store()

        if not np.all(adata.var_names == cell_state_df.index):
            raise ValueError(
                "adata.var_names should match cell_state_df.index, find interecting variables/genes first"
            )

        # add index for each cell (provided to pyro plate for correct minibatching)
        adata.obs["_indices"] = np.arange(adata.n_obs).astype("int64")
        scvi.data.register_tensor_from_anndata(
            adata,
            registry_key="ind_x",
            adata_attr_name="obs",
            adata_key_name="_indices",
        )

        super().__init__(adata)

        if model_class is None:
            model_class = LocationModelLinearDependentWMultiExperimentLocationBackgroundNormLevelGeneAlphaPyroModel

        self.cell_state_df_ = cell_state_df
        self.n_factors_ = cell_state_df.shape[1]
        self.factor_names_ = cell_state_df.columns.values

        if not detection_mean_per_sample:
            # compute expected change in sensitivity (m_g in V1 or y_s in V2)
            sc_total = cell_state_df.sum(0).mean()
            sp_total = get_from_registry(self.adata, _CONSTANTS.X_KEY).sum(1).mean()
            get_from_registry(adata, _CONSTANTS.BATCH_KEY)
            self.detection_mean_ = (
                sp_total / model_kwargs.get("N_cells_per_location", 1)
            ) / sc_total
            self.detection_mean_ = self.detection_mean_ * detection_mean_correction
            model_kwargs["detection_mean"] = self.detection_mean_
        else:
            # compute expected change in sensitivity (m_g in V1 and y_s in V2)
            sc_total = cell_state_df.sum(0).mean()
            sp_total = get_from_registry(self.adata, _CONSTANTS.X_KEY).sum(1)
            batch = get_from_registry(self.adata, _CONSTANTS.BATCH_KEY).flatten()
            sp_total = np.array(
                [
                    sp_total[batch == b].mean()
                    for b in range(self.summary_stats["n_batch"])
                ]
            )
            self.detection_mean_ = (
                sp_total / model_kwargs.get("N_cells_per_location", 1)
            ) / sc_total
            self.detection_mean_ = self.detection_mean_ * detection_mean_correction
            model_kwargs["detection_mean"] = self.detection_mean_.reshape(
                (self.summary_stats["n_batch"], 1)
            ).astype("float32")

        detection_alpha = model_kwargs.get("detection_alpha", None)
        if detection_alpha is not None:
            if type(detection_alpha) is dict:
                batch_mapping = self.adata.uns["_scvi"]["categorical_mappings"][
                    "_scvi_batch"
                ]["mapping"]
                self.detection_alpha_ = pd.Series(detection_alpha)[batch_mapping]
                model_kwargs["detection_alpha"] = self.detection_alpha_.values.reshape(
                    (self.summary_stats["n_batch"], 1)
                ).astype("float32")

        self.module = Cell2locationBaseModule(
            model=model_class,
            n_obs=self.summary_stats["n_cells"],
            n_vars=self.summary_stats["n_vars"],
            n_factors=self.n_factors_,
            n_batch=self.summary_stats["n_batch"],
            cell_state_mat=self.cell_state_df_.values.astype("float32"),
            **model_kwargs,
        )
        self._model_summary_string = f'cell2location model with the following params: \nn_factors: {self.n_factors_} \nn_batch: {self.summary_stats["n_batch"]} '
        self.init_params_ = self._get_init_params(locals())

    def export_posterior(
        self,
        adata,
        sample_kwargs: Optional[dict] = None,
        export_slot: str = "mod",
        add_to_obsm: list = ["means", "stds", "q05", "q95"],
    ):
        """
        Summarise posterior distribution and export results (cell abundance) to anndata object:
        1. adata.obsm: Estimated cell abundance as pd.DataFrames for each posterior distribution summary `add_to_obsm`,
            posterior mean, sd, 5% and 95% quantiles (['means', 'stds', 'q05', 'q95']).
            If export to adata.obsm fails with error, results are saved to adata.obs instead.
        2. adata.uns: Posterior of all parameters, model name, date,
            cell type names ('factor_names'), obs and var names.

        Parameters
        ----------
        adata
            anndata object where results should be saved
        sample_kwargs
            arguments for self.sample_posterior (generating and summarising posterior samples), namely:
                num_samples - number of samples to use (Default = 1000).
                batch_size - data batch size (keep low enough to fit on GPU, default 2048).
                use_gpu - use gpu for generating samples?
        export_slot
            adata.uns slot where to export results
        add_to_obsm
            posterior distribution summary to export in adata.obsm (['means', 'stds', 'q05', 'q95']).
        Returns
        -------

        """

        sample_kwargs = sample_kwargs if isinstance(sample_kwargs, dict) else dict()

        # generate samples from posterior distributions for all parameters
        # and compute mean, 5%/95% quantiles and standard deviation
        self.samples = self.sample_posterior(**sample_kwargs)

        # export posterior distribution summary for all parameters and
        # annotation (model, date, var, obs and cell type names) to anndata object
        adata.uns[export_slot] = self._export2adata(self.samples)

        # add estimated cell abundance as dataframe to obsm in anndata
        # first convert np.arrays to pd.DataFrames with cell type and observation names
        # data frames contain mean, 5%/95% quantiles and standard deviation, denoted by a prefix
        for k in add_to_obsm:
            sample_df = self.sample2df_obs(
                self.samples,
                site_name="w_sf",
                summary_name=k,
                name_prefix="cell_abundance",
            )
            try:
                adata.obsm[f"{k}_cell_abundance_w_sf"] = sample_df.loc[
                    adata.obs.index, :
                ]
            except ValueError:
                # Catching weird error with obsm: `ValueError: value.index does not match parent’s axis 1 names`
                adata.obs[sample_df.columns] = sample_df.loc[adata.obs.index, :]

        return adata
