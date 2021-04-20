import numpy as np
import pandas as pd
from anndata import AnnData
from scipy.sparse import csr_matrix

import scvi
from scvi.model.base import BaseModelClass

from ._base import Cell2locationTrainSampleMixin, PltExportMixin
from ._module import Cell2locationModule


def setup_anndata(adata, cell_state_df, **kwargs):
    """Setup anndata for reference-based decomposition:

    1) subset adata and cell_state_df to common variables
    2) add necesary information to adata.uns['_scvi']

    :param **kwargs: arguments for `scvi.data.setup_anndata`
    """

    intersect = np.intersect1d(cell_state_df.index, adata.var_names)
    cell_state_df = cell_state_df.loc[intersect, :]
    adata = adata[:, intersect]
    adata.varm["cell_state"] = cell_state_df.values.astype("float32")

    scvi.data.setup_anndata(adata, **kwargs)

    adata.uns["_scvi"]["summary_stats"]["n_obs"] = adata.n_obs

    # add cell type number and names
    adata.uns["_scvi"]["summary_stats"]["n_factors"] = cell_state_df.shape[1]
    adata.uns["_scvi"]["categorical_mappings"]["_scvi_factors"] = {
        "varm_key": "cell_state",
        "mapping": cell_state_df.columns.values,
    }

    # add index for each cell (provided to pyro plate for correct minibatching
    adata.obs["_indices"] = np.arange(adata.n_obs).astype("int64")
    adata.uns["_scvi"]["data_registry"]["ind_x"] = {
        "attr_name": "obs",
        "attr_key": "_indices",
    }

    return adata


def compute_cluster_averages(adata, labels, use_raw=True, layer=None):
    """
    :param adata: AnnData object of reference single-cell dataset
    :param labels: Name of adata.obs column containing cluster labels
    :returns: pd.DataFrame of cluster average expression of each gene
    """

    if layer is not None:
        x = adata.layers["layer"]
        var_names = adata.var_names
    else:
        if not use_raw:
            x = adata.X
            var_names = adata.var_names
        else:
            if not adata.raw:
                raise ValueError(
                    "AnnData object has no raw data, change `use_raw=True, layer=None` or fix your object"
                )
            x = adata.raw.X
            var_names = adata.raw.var_names

    if sum(adata.obs.columns == labels) != 1:
        raise ValueError("cluster_col is absent in adata_ref.obs or not unique")

    all_clusters = np.unique(adata.obs[labels])
    averages_mat = np.zeros((1, x.shape[1]))

    for c in all_clusters:
        sparse_subset = csr_matrix(x[np.isin(adata.obs[labels], c), :])
        aver = sparse_subset.mean(0)
        averages_mat = np.concatenate((averages_mat, aver))
    averages_mat = averages_mat[1:, :].T
    averages_df = pd.DataFrame(data=averages_mat, index=var_names, columns=all_clusters)

    return averages_df


class Cell2locationPltExportMixin:
    def _sample2df(self, node_name="w_sf"):
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
            self._sample2df()

        # add cell factors to adata
        adata.obs[self.w_sf_df.columns] = self.w_sf_df.loc[adata.obs.index, :]

        # add cell factor sd to adata
        adata.obs[self.w_sf_sd.columns] = self.w_sf_sd.loc[adata.obs.index, :]

        # add cell factor 5% and 95% quantiles to adata
        adata.obs[self.w_sf_q05.columns] = self.w_sf_q05.loc[adata.obs.index, :]
        adata.obs[self.w_sf_q95.columns] = self.w_sf_q95.loc[adata.obs.index, :]

        return adata


class Cell2location(
    Cell2locationTrainSampleMixin,
    BaseModelClass,
    PltExportMixin,
    Cell2locationPltExportMixin,
):
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
        batch_size=None,
        module=None,
        **model_kwargs,
    ):

        super().__init__(adata)

        if module is None:
            module = Cell2locationModule

        self.module = module(
            n_obs=self.summary_stats["n_obs"],
            n_vars=self.summary_stats["n_vars"],
            n_factors=self.summary_stats["n_factors"],
            n_batch=self.summary_stats["n_batch"],
            batch_size=batch_size,
            cell_state_mat=self.adata.varm[
                self.scvi_setup_dict_["categorical_mappings"]["_scvi_factors"][
                    "varm_key"
                ]
            ],
            **model_kwargs,
        )
