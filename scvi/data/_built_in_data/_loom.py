import logging
import os
import warnings

import numpy as np
import pandas as pd
from anndata import AnnData

from scvi.data._download import _download

logger = logging.getLogger(__name__)


def _load_retina(save_path: str = "data/") -> AnnData:
    """
    Loads retina dataset.

    The dataset of bipolar cells contains after their original pipeline for filtering 27,499 cells and
    13,166 genes coming from two batches. We use the cluster annotation from 15 cell-types from the author.
    We also extract their normalized data with Combat and use it for benchmarking.
    """
    save_path = os.path.abspath(save_path)
    url = "https://github.com/YosefLab/scVI-data/raw/master/retina.loom"
    save_fn = "retina.loom"
    _download(url, save_path, save_fn)
    adata = _load_loom(os.path.join(save_path, save_fn))
    cell_types = [
        "RBC",
        "MG",
        "BC5A",
        "BC7",
        "BC6",
        "BC5C",
        "BC1A",
        "BC3B",
        "BC1B",
        "BC2",
        "BC5D",
        "BC3A",
        "BC5B",
        "BC4",
        "BC8_9",
    ]
    adata.obs["labels"] = [
        cell_types[i] for i in adata.obs["ClusterID"].values.astype(int).ravel()
    ]
    del adata.obs["ClusterID"]
    adata.obs["batch"] = pd.Categorical(adata.obs["BatchID"].values.copy())
    del adata.obs["BatchID"]

    return adata


def _load_prefrontalcortex_starmap(save_path: str = "data/") -> AnnData:
    """Loads a starMAP dataset of 3,704 cells and 166 genes from the mouse pre-frontal cortex :cite:p:`Wang18`."""
    save_path = os.path.abspath(save_path)
    url = "https://github.com/YosefLab/scVI-data/raw/master/mpfc-starmap.loom"
    save_fn = "mpfc-starmap.loom"
    _download(url, save_path, save_fn)
    adata = _load_loom(os.path.join(save_path, save_fn))

    adata.obs["labels"] = adata.obs.Clusters.values
    del adata.obs["Clusters"]

    adata.obs["batch"] = adata.obs.BatchID.values
    del adata.obs["BatchID"]
    adata.obs["x_coord"] = adata.obsm["Spatial_coordinates"][:, 0]
    adata.obs["y_coord"] = adata.obsm["Spatial_coordinates"][:, 1]

    return adata


def _load_frontalcortex_dropseq(save_path: str = "data/") -> AnnData:
    save_path = os.path.abspath(save_path)
    url = "https://github.com/YosefLab/scVI-data/raw/master/fc-dropseq.loom"
    save_fn = "fc-dropseq.loom"
    _download(url, save_path, save_fn)
    adata = _load_loom(os.path.join(save_path, save_fn))
    adata.obs["batch"] = adata.obs["Clusters"]
    del adata.obs["Clusters"]
    adata.obs["labels"] = np.zeros(adata.shape[0], dtype=np.int64)

    # reorder labels such that layers of the cortex are in order
    # order_labels = [5, 6, 3, 2, 4, 0, 1, 8, 7, 9, 10, 11, 12, 13]
    # self.reorder_cell_types(self.cell_types[order_labels])

    return adata


def _load_annotation_simulation(name: str, save_path: str = "data/") -> AnnData:
    """
    Simulated datasets for scANVI tutorials.

    Parameters
    ----------
    name
        One of the following:
        * ``'1'``
        * ``'2'``
        * ``'3'``
    save_path
        Location for saving the dataset.
    """
    save_path = os.path.abspath(save_path)
    url = "https://github.com/YosefLab/scVI-data/raw/master/simulation/simulation_{}.loom".format(
        name
    )
    save_fn = f"simulation_{name}.loom"
    _download(url, save_path, save_fn)
    adata = _load_loom(os.path.join(save_path, save_fn))

    adata.obs["labels"] = adata.obs.ClusterID.values
    del adata.obs["ClusterID"]

    adata.obs["batch"] = adata.obs.BatchID.values
    del adata.obs["BatchID"]

    return adata


def _load_loom(path_to_file: str, gene_names_attribute_name: str = "Gene") -> AnnData:
    import loompy

    dataset = loompy.connect(path_to_file)
    select = dataset[:, :].sum(axis=0) > 0  # Take out cells that don't express any gene
    if not all(select):
        warnings.warn("Removing empty cells")

    var_dict, obs_dict, uns_dict, obsm_dict = {}, {}, {}, {}
    for row_key in dataset.ra:
        if row_key == gene_names_attribute_name:
            gene_names = dataset.ra[gene_names_attribute_name].astype(str)
        else:
            var_dict[row_key] = dataset.ra[row_key]
            if type(var_dict[row_key]) is np.ndarray:
                var_dict[row_key] = var_dict[row_key].ravel()

    for column_key in dataset.ca:
        obs_dict = obs_dict if obs_dict is not None else {}
        obs_dict[column_key] = dataset.ca[column_key][select]
        if type(obs_dict[column_key]) is np.ndarray:
            if len(obs_dict[column_key]) == len(obs_dict[column_key].ravel()):
                obs_dict[column_key] = obs_dict[column_key].ravel()
            else:
                obsm_dict[column_key] = obs_dict[column_key]
                del obs_dict[column_key]

    for global_key in dataset.attrs:
        uns_dict = uns_dict if uns_dict is not None else {}
        uns_dict[global_key] = dataset.attrs[global_key]
        if type(uns_dict[global_key]) is np.ndarray:
            uns_dict[global_key] = uns_dict[global_key].ravel()
    data = dataset[:, :].T  # change matrix to cells by genes
    dataset.close()

    adata = AnnData(X=data, obs=obs_dict, var=var_dict, uns=uns_dict, obsm=obsm_dict)
    adata = adata[select].copy()
    adata.var_names = gene_names

    return adata
