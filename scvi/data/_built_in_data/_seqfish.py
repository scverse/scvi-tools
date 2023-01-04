import logging
import os
import zipfile

import anndata
import numpy as np
import pandas as pd

from scvi.data._download import _download

logger = logging.getLogger(__name__)


def _load_seqfishplus(
    save_path: str = "data/",
    tissue_region: str = "subventricular cortex",
) -> anndata.AnnData:

    if tissue_region == "subventricular cortex":
        file_prefix = "cortex_svz"
    elif tissue_region == "olfactory bulb":
        file_prefix = "ob"
    else:
        raise ValueError(
            '`tissue_type` must be "subventricular cortex" or "olfactory bulb", but got {}'.format(
                tissue_region
            )
        )

    save_path = os.path.abspath(save_path)
    url = "https://github.com/CaiGroup/seqFISH-PLUS/raw/master/sourcedata.zip"
    save_fn = "seqfishplus.zip"

    _download(url, save_path, save_fn)
    adata = _load_seqfishplus_data(
        os.path.join(save_path, save_fn), file_prefix, save_path, gene_by_cell=False
    )
    adata.obs["batch"] = np.zeros(adata.shape[0], dtype=np.int64)
    adata.obs["labels"] = np.zeros(adata.shape[0], dtype=np.int64)

    return adata


def _load_seqfishplus_data(
    path_to_file: str, file_prefix: str, save_path: str, gene_by_cell: bool = False
) -> anndata.AnnData:
    counts_filename = f"sourcedata/{file_prefix}_counts.csv"
    coordinates_filename = f"sourcedata/{file_prefix}_cellcentroids.csv"
    extract_location = os.path.join(save_path, "seqfishplus")
    if not os.path.exists(extract_location):
        os.makedirs(extract_location)
    with zipfile.ZipFile(path_to_file) as f:
        f.extract(counts_filename, path=extract_location)
        f.extract(coordinates_filename, path=extract_location)

    df_counts = pd.read_csv(os.path.join(extract_location, counts_filename))
    adata = anndata.AnnData(df_counts)
    adata.var_names = df_counts.columns
    df_coordinates = pd.read_csv(os.path.join(extract_location, coordinates_filename))

    adata.obs["X"] = df_coordinates["X"].values
    adata.obs["Y"] = df_coordinates["Y"].values
    adata.obs["cell_id"] = df_coordinates["Cell ID"].values
    adata.obs["field_of_view"] = df_coordinates["Field of View"].values

    return adata


def _load_seqfish(save_path: str = "data/") -> anndata.AnnData:
    save_path = os.path.abspath(save_path)
    url = "https://www.cell.com/cms/attachment/2080562255/2072099886/mmc6.xlsx"
    save_fn = "SeqFISH.xlsx"
    _download(url, save_path, save_fn)
    adata = _load_seqfish_data(os.path.join(save_path, save_fn))
    adata.obs["batch"] = np.zeros(adata.shape[0], dtype=np.int64)
    adata.obs["labels"] = np.zeros(adata.shape[0], dtype=np.int64)
    return adata


def _load_seqfish_data(path_to_file: str) -> anndata.AnnData:
    logger.info(f"Loading seqfish dataset from {path_to_file}")
    counts = pd.read_excel(
        path_to_file, sheet_name="Hippocampus Counts", engine="openpyxl"
    )
    data = (
        counts.values[:, 1:].astype(int).T
    )  # transpose because counts is genes X cells
    gene_names = counts.values[:, 0].astype(str)
    adata = anndata.AnnData(pd.DataFrame(data=data, columns=gene_names))
    logger.info("Finished loading seqfish dataset")
    return adata
