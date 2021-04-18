import numpy as np
import pandas as pd
from _base import PltExportMixin, TrainSampleMixin
from anndata import AnnData

# from scvi.data import register_tensor_from_anndata
from scvi.external.cell2location._module import Cell2locationModule
from scvi.model.base import BaseModelClass


class Cell2locationPltExportMixin:
    def sample2df(self, node_name="w_sf"):
        r"""Export spatial cell abundance as Pandas data frames.

        :param node_name: name of the model parameter to be exported
        :return: 4 Pandas dataframes added to model object:
            .spot_factors_df, .spot_factors_sd, .spot_factors_q05, .spot_factors_q95
        """

        if len(self.samples) == 0:
            raise ValueError(
                "Please run `.sample_posterior()` first to generate samples & summarise posterior of each parameter"
            )

        self.w_sf_df = pd.DataFrame.from_records(
            self.samples["post_sample_means"][node_name],
            index=self.obs_names,
            columns=["mean_" + node_name + i for i in self.fact_names],
        )

        self.w_sf_sd = pd.DataFrame.from_records(
            self.samples["post_sample_sds"][node_name],
            index=self.obs_names,
            columns=["sd_" + node_name + i for i in self.fact_names],
        )

        self.w_sf_q05 = pd.DataFrame.from_records(
            self.samples["post_sample_q05"][node_name],
            index=self.obs_names,
            columns=["q05_" + node_name + i for i in self.fact_names],
        )

        self.w_sf_q95 = pd.DataFrame.from_records(
            self.samples["post_sample_q95"][node_name],
            index=self.obs_names,
            columns=["q95_" + node_name + i for i in self.fact_names],
        )

    def annotate_spatial_adata(self, adata):
        r"""Add spatial cell abundance to adata.obs

        :param adata: anndata object to annotate
        :return: updated anndata object
        """

        if self.spot_factors_df is None:
            self.sample2df()

        # add cell factors to adata
        adata.obs[self.w_sf_df.columns] = self.w_sf_df.loc[adata.obs.index, :]

        # add cell factor sd to adata
        adata.obs[self.w_sf_sd.columns] = self.w_sf_sd.loc[adata.obs.index, :]

        # add cell factor 5% and 95% quantiles to adata
        adata.obs[self.w_sf_q05.columns] = self.w_sf_q05.loc[adata.obs.index, :]
        adata.obs[self.w_sf_q95.columns] = self.w_sf_q95.loc[adata.obs.index, :]

        return adata


class Cell2locationBaseModelClass(
    BaseModelClass, TrainSampleMixin, PltExportMixin, Cell2locationPltExportMixin
):
    """
    Reimplementation of cell2location [Kleshchevnikov20]_ model. Cell2locationBaseModelClass.

    https://github.com/BayraktarLab/cell2location

    Parameters
    ----------
    sc_adata
        single-cell AnnData object that has been registered via :func:`~scvi.data.setup_anndata`.
    use_gpu
        Use the GPU or not.
    **model_kwargs
        Keyword args for :class:`~scvi.external.cell2location...`

    Examples
    --------
    >>>
    """

    def __init__(
        self,
        adata: AnnData,
        cell_state_df: pd.DataFrame,
        var_names_read=None,
        sample_id=None,
        use_gpu: bool = True,
        batch_size: int = 1024,
        **model_kwargs,
    ):

        intersect = np.intersect1d(cell_state_df.index, adata.var_names)
        cell_state_df = cell_state_df.loc[intersect, :]
        adata = adata[:, intersect]
        adata.varm["cell_state_varm"] = cell_state_df.values

        super().__init__(
            adata,
            n_fact=cell_state_df.shape[1],
            var_names_read=var_names_read,
            fact_names=cell_state_df.columns,
            sample_id=sample_id,
            use_gpu=use_gpu,
            batch_size=batch_size,
        )

        self.cell_state_df = cell_state_df


class Cell2location(Cell2locationBaseModelClass):
    """
    Reimplementation of cell2location [Kleshchevnikov20]_ model. User-end model class.

    https://github.com/BayraktarLab/cell2location

    Parameters
    ----------
    sc_adata
        single-cell AnnData object that has been registered via :func:`~scvi.data.setup_anndata`.
    use_gpu
        Use the GPU or not.
    **model_kwargs
        Keyword args for :class:`~scvi.external.cell2location...`

    Examples
    --------
    >>>
    """

    def __init__(
        self,
        adata: AnnData,
        cell_state_df: pd.DataFrame,
        var_names_read=None,
        sample_id=None,
        module=None,
        use_gpu: bool = True,
        batch_size: int = 2048,
        **model_kwargs,
    ):

        super(Cell2location, self).__init__(
            adata=adata,
            cell_state_df=cell_state_df,
            var_names_read=var_names_read,
            sample_id=sample_id,
            use_gpu=use_gpu,
            batch_size=batch_size,
        )

        if module is None:
            module = Cell2locationModule

        self.scvi_setup_dict_

        self.module = module(
            n_obs=self.n_obs,
            n_var=self.n_var,
            n_fact=self.n_fact,
            n_exper=self.n_exper,
            batch_size=self.batch_size,
            cell_state_mat=self.cell_state_df.values,
            **model_kwargs,
        )
