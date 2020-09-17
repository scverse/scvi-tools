import pandas as pd
import numpy as np
import anndata
import os

from scvi.data import setup_anndata
from scvi.data._built_in_data._utils import _download


def _load_pbmcs_10x_cite_seq(
    save_path: str = "data/",
    protein_join: str = "inner",
    run_setup_anndata: bool = True,
):
    """Filtered PBMCs from 10x Genomics profiled with RNA and protein

    Datasets were filtered for doublets and other outliers as in
    https://github.com/YosefLab/totalVI_reproducibility/blob/master/data/data_filtering_scripts/pbmc_10k/pbmc_10k.py

    Parameters
    ----------
    protein_join
        Whether to take an inner join or outer join of proteins

    Returns
    -------
    `AnnData` with `.obsm["protein_expression"]

    Missing protein values are zero, and are identified during `AnnData` setup.
    """

    url = "https://github.com/YosefLab/scVI-data/raw/master/pbmc_10k_protein_v3.h5ad?raw=true"
    save_fn = "pbmc_10k_protein_v3.h5ad"
    _download(url, save_path, save_fn)
    dataset1 = anndata.read_h5ad(os.path.join(save_path, save_fn))

    url = "https://github.com/YosefLab/scVI-data/raw/master/pbmc_5k_protein_v3.h5ad?raw=true"
    save_fn = "pbmc_5k_protein_v3.h5ad"
    _download(url, save_path, save_fn)
    dataset2 = anndata.read_h5ad(os.path.join(save_path, "pbmc_5k_protein_v3.h5ad"))

    common_genes = dataset1.var_names.intersection(dataset2.var_names)
    dataset1 = dataset1[:, common_genes]
    dataset2 = dataset2[:, common_genes]
    dataset1.obsm["protein_expression"] = pd.DataFrame(
        dataset1.obsm["protein_expression"],
        columns=dataset1.uns["protein_names"],
        index=dataset1.obs_names,
    )
    dataset2.obsm["protein_expression"] = pd.DataFrame(
        dataset2.obsm["protein_expression"],
        columns=dataset2.uns["protein_names"],
        index=dataset2.obs_names,
    )
    del dataset1.uns["protein_names"]
    del dataset2.uns["protein_names"]

    dataset = dataset1.concatenate(dataset2, join=protein_join)
    dataset.obsm["protein_expression"] = dataset.obsm["protein_expression"].fillna(0)
    dataset.obs["labels"] = np.zeros(dataset.shape[0], dtype=np.int64)
    dataset.obs["batch"] = dataset.obs["batch"].astype(np.int64)

    if run_setup_anndata:
        setup_anndata(
            dataset,
            batch_key="batch",
            labels_key="labels",
            protein_expression_obsm_key="protein_expression",
        )

    return dataset
