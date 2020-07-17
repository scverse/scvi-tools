import logging
import os
import pandas as pd
import loompy
import numpy as np
import anndata

from scvi.dataset._utils import _download
from scvi.dataset import setup_anndata

logger = logging.getLogger(__name__)

_subtype_to_high_level_mapping = {
    "Astrocytes": ("Astrocyte Gfap", "Astrocyte Mfge8"),
    "Endothelials": ("Endothelial", "Endothelial 1"),
    "Inhibitory": (
        "Inhibitory Cnr1",
        "Inhibitory Kcnip2",
        "Inhibitory Pthlh",
        "Inhibitory Crhbp",
        "Inhibitory CP",
        "Inhibitory IC",
        "Inhibitory Vip",
    ),
    "Microglias": ("Perivascular Macrophages", "Microglia"),
    "Oligodendrocytes": (
        "Oligodendrocyte Precursor cells",
        "Oligodendrocyte COP",
        "Oligodendrocyte NF",
        "Oligodendrocyte Mature",
        "Oligodendrocyte MF",
    ),
    "Pyramidals": (
        "Pyramidal L2-3",
        "Pyramidal Cpne5",
        "Pyramidal L2-3 L5",
        "pyramidal L4",
        "Pyramidal L3-4",
        "Pyramidal Kcnip2",
        "Pyramidal L6",
        "Pyramidal L5",
        "Hippocampus",
    ),
}


def smfish(save_path="data/", use_high_level_cluster=True):
    save_path = os.path.abspath(save_path)
    url = "http://linnarssonlab.org/osmFISH/osmFISH_SScortex_mouse_all_cells.loom"
    save_fn = "osmFISH_SScortex_mouse_all_cell.loom"
    _download(url, save_path, save_fn)
    adata = _load_smfish(
        os.path.join(save_path, save_fn), use_high_level_cluster=use_high_level_cluster
    )
    return adata


def _load_smfish(path_to_file, use_high_level_cluster):
    logger.info("Loading smFISH dataset")
    ds = loompy.connect(path_to_file)
    x_coord, y_coord = ds.ca["X"], ds.ca["Y"]
    data = ds[:, :].T
    gene_names = ds.ra["Gene"].astype(np.str)
    labels = ds.ca["ClusterID"]
    str_labels = np.asarray(ds.ca["ClusterName"])
    labels_mapping = pd.Categorical(str_labels).categories

    if use_high_level_cluster:
        for high_level_cluster, subtypes in _subtype_to_high_level_mapping.items():
            for subtype in subtypes:
                idx = np.where(str_labels == subtype)
                str_labels[idx] = high_level_cluster
        cell_types_to_keep = [
            "Astrocytes",
            "Endothelials",
            "Inhibitory",
            "Microglias",
            "Oligodendrocytes",
            "Pyramidals",
        ]
        row_indices = [
            i
            for i in range(data.shape[0])
            if ds.ca["ClusterName"][i] in cell_types_to_keep
        ]
        str_labels = str_labels[row_indices]
        data = data[row_indices, :]
        x_coord = x_coord[row_indices]
        y_coord = y_coord[row_indices]

        str_labels = pd.Categorical(str_labels)
        labels = str_labels.codes
        labels_mapping = str_labels.categories

    adata = anndata.AnnData(
        X=data,
        obs={
            "x_coord": x_coord,
            "y_coord": y_coord,
            "labels": labels,
            "str_labels": str_labels,
        },
        uns={"cell_types": labels_mapping},
    )
    adata.var_names = gene_names
    return adata
