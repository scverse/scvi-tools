import csv
import logging
import os

import anndata
import numpy as np
import pandas as pd

from scvi.data._anndata import setup_anndata
from scvi.data._built_in_data._download import _download

logger = logging.getLogger(__name__)


def _load_cortex(
    save_path: str = "data/", run_setup_anndata: bool = True
) -> anndata.AnnData:
    """Loads cortex dataset."""
    save_path = os.path.abspath(save_path)
    url = "https://storage.googleapis.com/linnarsson-lab-www-blobs/blobs/cortex/expression_mRNA_17-Aug-2014.txt"
    save_fn = "expression.bin"
    _download(url, save_path, save_fn)
    adata = _load_cortex_txt(os.path.join(save_path, save_fn))
    if run_setup_anndata:
        setup_anndata(adata, labels_key="labels")
    return adata


def _load_cortex_txt(path_to_file: str) -> anndata.AnnData:
    logger.info("Loading Cortex data from {}".format(path_to_file))
    rows = []
    gene_names = []
    with open(path_to_file, "r") as csvfile:
        data_reader = csv.reader(csvfile, delimiter="\t")
        for i, row in enumerate(data_reader):
            if i == 1:
                precise_clusters = np.asarray(row, dtype=str)[2:]
            if i == 8:
                clusters = np.asarray(row, dtype=str)[2:]
            if i >= 11:
                rows.append(row[1:])
                gene_names.append(row[0])
    cell_types, labels = np.unique(clusters, return_inverse=True)
    _, precise_labels = np.unique(precise_clusters, return_inverse=True)
    data = np.asarray(rows, dtype=np.int).T[1:]
    gene_names = np.asarray(gene_names, dtype=np.str)
    gene_indices = []

    extra_gene_indices = []
    gene_indices = np.concatenate([gene_indices, extra_gene_indices]).astype(np.int32)
    if gene_indices.size == 0:
        gene_indices = slice(None)

    data = data[:, gene_indices]
    gene_names = gene_names[gene_indices]
    data_df = pd.DataFrame(data, columns=gene_names)
    adata = anndata.AnnData(X=data_df)
    adata.obs["labels"] = labels
    adata.obs["precise_labels"] = precise_clusters
    adata.obs["cell_type"] = clusters
    logger.info("Finished loading Cortex data")
    return adata
