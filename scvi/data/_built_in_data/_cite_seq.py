import os

import anndata
import pandas as pd

from scvi.data import setup_anndata
from scvi.data._built_in_data._download import _download


def _load_pbmcs_10x_cite_seq(
    save_path: str = "data/",
    protein_join: str = "inner",
    run_setup_anndata: bool = True,
):
    """
    Filtered PBMCs from 10x Genomics profiled with RNA and protein.

    Datasets were filtered for doublets and other outliers as in
    https://github.com/YosefLab/totalVI_reproducibility/blob/master/data/data_filtering_scripts/pbmc_10k/pbmc_10k.py

    Parameters
    ----------
    save_path
        Location to use when saving/loading the data.
    protein_join
        Whether to take an inner join or outer join of proteins
    run_setup_anndata
        If true, runs setup_anndata() on dataset before returning

    Returns
    -------
    `AnnData` with `.obsm["protein_expression"]

    Missing protein values are zero, and are identified during `AnnData` setup.
    """
    url = "https://github.com/YosefLab/scVI-data/raw/master/pbmc_10k_protein_v3.h5ad?raw=true"
    save_fn = "pbmc_10k_protein_v3.h5ad"
    _download(url, save_path, save_fn)
    dataset1 = anndata.read_h5ad(os.path.join(save_path, save_fn))
    dataset1.obs["batch"] = "PBMC10k"

    url = "https://github.com/YosefLab/scVI-data/raw/master/pbmc_5k_protein_v3.h5ad?raw=true"
    save_fn = "pbmc_5k_protein_v3.h5ad"
    _download(url, save_path, save_fn)
    dataset2 = anndata.read_h5ad(os.path.join(save_path, "pbmc_5k_protein_v3.h5ad"))
    dataset2.obs["batch"] = "PBMC5k"

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

    dataset = anndata.concat([dataset1, dataset2], join=protein_join)
    dataset.obsm["protein_expression"] = dataset.obsm["protein_expression"].fillna(0)

    if run_setup_anndata:
        setup_anndata(
            dataset,
            batch_key="batch",
            protein_expression_obsm_key="protein_expression",
        )

    return dataset


def _load_spleen_lymph_cite_seq(
    save_path: str = "data/",
    protein_join: str = "inner",
    remove_outliers: bool = True,
    run_setup_anndata: bool = True,
):
    """
    Immune cells from the murine spleen and lymph nodes [GayosoSteier20]_.

    This dataset was used throughout the totalVI manuscript, and named SLN-all.

    Parameters
    ----------
    save_path
        Location to use when saving/loading the data.
    protein_join
        Whether to take an inner join or outer join of proteins
    remove_outliers
        Whether to remove clusters annotated as doublet or low quality
    run_setup_anndata
        If true, runs setup_anndata() on dataset before returning

    Returns
    -------
    `AnnData` with `.obsm["protein_expression"]

    Missing protein values are zero, and are identified during `AnnData` setup.
    """
    url = "https://github.com/YosefLab/scVI-data/raw/master/sln_111.h5ad?raw=true"
    save_fn = "sln_111.h5ad"
    _download(url, save_path, save_fn)
    dataset1 = anndata.read_h5ad(os.path.join(save_path, save_fn))
    dataset1.obsm["isotypes_htos"] = dataset1.obsm["htos"].copy()
    del dataset1.obsm["htos"]

    url = "https://github.com/YosefLab/scVI-data/raw/master/sln_208.h5ad?raw=true"
    save_fn = "sln_208.h5ad"
    _download(url, save_path, save_fn)
    dataset2 = anndata.read_h5ad(os.path.join(save_path, save_fn))

    common_genes = dataset1.var_names.intersection(dataset2.var_names)
    dataset1 = dataset1[:, common_genes]
    dataset2 = dataset2[:, common_genes]

    del dataset1.uns["protein_names"]
    del dataset2.uns["protein_names"]

    dataset = anndata.concat(
        [dataset1, dataset2],
        join=protein_join,
    )
    dataset.obsm["protein_expression"] = dataset.obsm["protein_expression"].fillna(0)

    if remove_outliers:
        include_cells = [
            c not in ["16,0", "17", "19", "21", "23", "24,0", "24,2", "25", "29"]
            for c in dataset.obs["leiden_subclusters"]
        ]
        dataset = dataset[include_cells].copy()

    if run_setup_anndata:
        setup_anndata(
            dataset,
            batch_key="batch",
            labels_key="cell_types",
            protein_expression_obsm_key="protein_expression",
        )

    return dataset
