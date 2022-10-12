import logging
import os

import anndata
import numpy as np

from scvi.data._download import _download

logger = logging.getLogger(__name__)


def _load_breast_cancer_dataset(save_path: str = "data/"):
    save_path = os.path.abspath(save_path)
    url = "https://www.spatialresearch.org/wp-content/uploads/2016/07/Layer2_BC_count_matrix-1.tsv"
    save_fn = "Layer2_BC_count_matrix-1.tsv"
    _download(url, save_path, save_fn)
    adata = _load_csv(
        os.path.join(save_path, save_fn), delimiter="\t", gene_by_cell=False
    )
    adata.obs["batch"] = np.zeros(adata.shape[0]).astype(int)
    adata.obs["labels"] = np.zeros(adata.shape[0]).astype(int)

    return adata


def _load_mouse_ob_dataset(save_path: str = "data/"):
    save_path = os.path.abspath(save_path)
    url = "https://www.spatialresearch.org/wp-content/uploads/2016/07/Rep11_MOB_count_matrix-1.tsv"
    save_fn = "Rep11_MOB_count_matrix-1.tsv"
    _download(url, save_path, save_fn)
    adata = _load_csv(
        os.path.join(save_path, save_fn), delimiter="\t", gene_by_cell=False
    )
    adata.obs["batch"] = np.zeros(adata.shape[0]).astype(int)
    adata.obs["labels"] = np.zeros(adata.shape[0]).astype(int)

    return adata


def _load_csv(
    path_to_file: str,
    gene_by_cell: bool = False,
    delimiter: str = ",",
    first_column_names: bool = None,
):
    logger.info(f"Loading dataset from {path_to_file}")
    adata = anndata.read_csv(
        path_to_file, delimiter=delimiter, first_column_names=first_column_names
    )
    if gene_by_cell:
        adata.X = adata.X.T
    logger.info("Finished loading dataset")
    return adata
