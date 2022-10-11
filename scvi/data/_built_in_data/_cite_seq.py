import os

import anndata
import numpy as np
import pandas as pd

from scvi import settings
from scvi.data._download import _download


def _load_pbmcs_10x_cite_seq(
    save_path: str = "data/",
    protein_join: str = "inner",
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

    Returns
    -------
    ``AnnData`` with ``.obsm['protein_expression']``
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

    return dataset


def _load_spleen_lymph_cite_seq(
    save_path: str = "data/",
    protein_join: str = "inner",
    remove_outliers: bool = True,
):
    """
    Immune cells from the murine spleen and lymph nodes :cite:p:`GayosoSteier21`.

    This dataset was used throughout the totalVI manuscript, and named SLN-all.

    Parameters
    ----------
    save_path
        Location to use when saving/loading the data.
    protein_join
        Whether to take an inner join or outer join of proteins
    remove_outliers
        Whether to remove clusters annotated as doublet or low quality

    Returns
    -------
    ``AnnData`` with ``.obsm['protein_expression']``
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

    return dataset


def _load_pbmc_seurat_v4_cite_seq(
    save_path: str = "data/",
    apply_filters: bool = True,
    aggregate_proteins: bool = True,
    mask_protein_batches: int = 0,
):
    url = "https://ndownloader.figshare.com/files/27458840"
    save_fn = "pbmc_seurat_v4.h5ad"
    _download(url, save_path, save_fn)
    adata = anndata.read_h5ad(os.path.join(save_path, save_fn))

    if aggregate_proteins:
        protein_dict = {}
        ref_proteins = adata.obsm["protein_counts"].columns
        for p in ref_proteins:
            if p.split("-")[-1] == "1" or p.split("-")[-1] == "2":
                root = p.split("-")[0]
                if root not in ["Notch", "TCR"]:
                    try:
                        protein_dict[root] = np.asarray(
                            adata.obsm["protein_counts"][root + "-1"]
                            + adata.obsm["protein_counts"][root + "-2"]
                        )
                    except KeyError:
                        protein_dict[p] = np.asarray(adata.obsm["protein_counts"][p])
                else:
                    protein_dict[p] = np.asarray(adata.obsm["protein_counts"][p])
            else:
                protein_dict[p] = np.asarray(adata.obsm["protein_counts"][p])
        protein_df = pd.DataFrame.from_dict(protein_dict)
        protein_df.index = adata.obsm["protein_counts"].index
        adata.obsm["protein_counts"] = protein_df

    if apply_filters:
        adata.obs["total_counts"] = np.ravel(adata.X.sum(axis=1).A)
        adata.var["mt"] = adata.var_names.str.startswith("MT-")
        adata.obs["total_counts_mt"] = np.ravel(
            adata.X[:, adata.var["mt"].values].sum(axis=1).A
        )
        adata.obs["pct_counts_mt"] = (
            adata.obs["total_counts_mt"] / adata.obs["total_counts"] * 100
        )

        adata.obs["Protein log library size"] = np.log(
            adata.obsm["protein_counts"].sum(1)
        )
        adata.obs["Number proteins detected"] = (adata.obsm["protein_counts"] > 0).sum(
            1
        )
        adata.obs["RNA log library size"] = np.log(adata.X.sum(1).A)

        # actually filter
        adata = adata[adata.obs["Protein log library size"] > 7.6]
        adata = adata[adata.obs["Protein log library size"] < 10.3]
        adata = adata[adata.obs["Number proteins detected"] > 150]
        # filter doublet
        adata = adata[adata.obs["celltype.l2"] != "Doublet"]
        # MT
        adata = adata[adata.obs["pct_counts_mt"] < 12].copy()

    if mask_protein_batches > 24:
        raise ValueError("mask_protein_batches must be less than 24")

    if mask_protein_batches > 0:
        random_state = np.random.RandomState(seed=settings.seed)
        rand_cats = random_state.permutation(
            adata.obs["orig.ident"].astype("category").cat.categories
        )[:mask_protein_batches]
        for r in rand_cats:
            adata.obsm["protein_counts"][adata.obs["orig.ident"] == r] = 0.0

    return adata
