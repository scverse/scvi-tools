import logging
import os

import numpy as np
import pandas as pd
from anndata import AnnData
from anndata.io import read_loom

from scvi.data._download import _download

logger = logging.getLogger(__name__)


def _load_retina(save_path: str = "data/") -> AnnData:
    """Loads retina dataset.

    The dataset of bipolar cells contains after their original pipeline for filtering 27,499 cells
    and 13,166 genes coming from two batches. We use the cluster annotation from 15 cell-types from
    the author. We also extract their normalized data with Combat and use it for benchmarking.
    """
    save_path = os.path.abspath(save_path)
    url = "https://github.com/YosefLab/scVI-data/raw/master/retina.loom"
    save_fn = "retina.loom"
    _download(url, save_path, save_fn)
    # Check numpy version for loompy
    if np.__version__ >= "2.0.0":
        raise ValueError("reading a loom file requires Numpy version smaller than 2.0.0")
    adata = read_loom(os.path.join(save_path, save_fn))
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
    """Loads a starMAP dataset from the mouse pre-frontal cortex :cite:p:`Wang18`.

    Contains 3,704 cells and 166 genes.
    """
    save_path = os.path.abspath(save_path)
    url = "https://github.com/YosefLab/scVI-data/raw/master/mpfc-starmap.loom"
    save_fn = "mpfc-starmap.loom"
    _download(url, save_path, save_fn)
    # Check numpy version for loompy
    if np.__version__ >= "2.0.0":
        raise ValueError("reading a loom file requires Numpy version smaller than 2.0.0")
    adata = read_loom(os.path.join(save_path, save_fn))

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
    # Check numpy version for loompy
    if np.__version__ >= "2.0.0":
        raise ValueError("reading a loom file requires Numpy version smaller than 2.0.0")
    adata = read_loom(os.path.join(save_path, save_fn))
    adata.obs["batch"] = adata.obs["Clusters"]
    del adata.obs["Clusters"]
    adata.obs["labels"] = np.zeros(adata.shape[0], dtype=np.int64)

    # reorder labels such that layers of the cortex are in order
    # order_labels = [5, 6, 3, 2, 4, 0, 1, 8, 7, 9, 10, 11, 12, 13]
    # self.reorder_cell_types(self.cell_types[order_labels])

    return adata


def _load_annotation_simulation(name: str, save_path: str = "data/") -> AnnData:
    """Simulated datasets for scANVI tutorials.

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
    url = f"https://github.com/YosefLab/scVI-data/raw/master/simulation/simulation_{name}.loom"
    save_fn = f"simulation_{name}.loom"
    _download(url, save_path, save_fn)
    if np.__version__ >= "2.0.0":
        raise ValueError("reading a loom file requires Numpy version smaller than 2.0.0")
    adata = read_loom(os.path.join(save_path, save_fn))

    adata.obs["labels"] = adata.obs.ClusterID.values
    del adata.obs["ClusterID"]

    adata.obs["batch"] = adata.obs.BatchID.values
    del adata.obs["BatchID"]

    return adata
