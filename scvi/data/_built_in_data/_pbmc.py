import os
import pickle
from typing import List

import anndata
import numpy as np
import pandas as pd

from scvi.data import setup_anndata
from scvi.data._built_in_data._dataset_10x import _load_dataset_10x
from scvi.data._built_in_data._download import _download


def _load_purified_pbmc_dataset(
    save_path: str = "data/",
    subset_datasets: List[str] = None,
    run_setup_anndata: bool = True,
) -> anndata.AnnData:
    url = "https://github.com/YosefLab/scVI-data/raw/master/PurifiedPBMCDataset.h5ad"
    save_fn = "PurifiedPBMCDataset.h5ad"
    _download(url, save_path, save_fn)
    path_to_file = os.path.join(save_path, save_fn)
    adata = anndata.read(path_to_file)

    dataset_names = [
        "cd4_t_helper",
        "regulatory_t",
        "naive_t",
        "memory_t",
        "cytotoxic_t",
        "naive_cytotoxic",
        "b_cells",
        "cd4_t_helper",
        "cd34",
        "cd56_nk",
        "cd14_monocytes",
    ]
    if subset_datasets is not None:
        row_indices = []
        for dataset in subset_datasets:
            assert dataset in dataset_names
            idx = np.where(adata.obs["cell_types"] == dataset)[0]
            row_indices.append(idx)
        row_indices = np.concatenate(row_indices)
        adata = adata[row_indices].copy()

    if run_setup_anndata:
        setup_anndata(adata, batch_key="batch", labels_key="labels")

    return adata


def _load_pbmc_dataset(
    save_path: str = "data/",
    run_setup_anndata: bool = True,
    remove_extracted_data: bool = True,
) -> anndata.AnnData:
    urls = [
        "https://github.com/YosefLab/scVI-data/raw/master/gene_info.csv",
        "https://github.com/YosefLab/scVI-data/raw/master/pbmc_metadata.pickle",
    ]
    save_fns = ["gene_info_pbmc.csv", "pbmc_metadata.pickle"]

    for i in range(len(urls)):
        _download(urls[i], save_path, save_fns[i])

    de_metadata = pd.read_csv(os.path.join(save_path, "gene_info_pbmc.csv"), sep=",")
    pbmc_metadata = pickle.load(
        open(os.path.join(save_path, "pbmc_metadata.pickle"), "rb")
    )
    pbmc8k = _load_dataset_10x(
        "pbmc8k",
        save_path=save_path,
        var_names="gene_ids",
        remove_extracted_data=remove_extracted_data,
    )
    pbmc4k = _load_dataset_10x(
        "pbmc4k",
        save_path=save_path,
        var_names="gene_ids",
        remove_extracted_data=remove_extracted_data,
    )
    barcodes = np.concatenate((pbmc8k.obs_names, pbmc4k.obs_names))

    adata = pbmc8k.concatenate(pbmc4k)
    adata.obs_names = barcodes

    dict_barcodes = dict(zip(barcodes, np.arange(len(barcodes))))
    subset_cells = []
    barcodes_metadata = pbmc_metadata["barcodes"].index.values.ravel().astype(np.str)
    for barcode in barcodes_metadata:
        if (
            barcode in dict_barcodes
        ):  # barcodes with end -11 filtered on 10X website (49 cells)
            subset_cells += [dict_barcodes[barcode]]
    adata = adata[np.asarray(subset_cells), :].copy()
    idx_metadata = np.asarray(
        [not barcode.endswith("11") for barcode in barcodes_metadata], dtype=np.bool
    )
    genes_to_keep = list(
        de_metadata["ENSG"].values
    )  # only keep the genes for which we have de data
    difference = list(
        set(genes_to_keep).difference(set(adata.var_names))
    )  # Non empty only for unit tests
    for gene in difference:
        genes_to_keep.remove(gene)

    adata = adata[:, genes_to_keep].copy()
    design = pbmc_metadata["design"][idx_metadata]
    raw_qc = pbmc_metadata["raw_qc"][idx_metadata]
    normalized_qc = pbmc_metadata["normalized_qc"][idx_metadata]

    design.index = adata.obs_names
    raw_qc.index = adata.obs_names
    normalized_qc.index = adata.obs_names
    adata.obs["batch"] = adata.obs["batch"].astype(np.int64)
    adata.obsm["design"] = design
    adata.obsm["raw_qc"] = raw_qc
    adata.obsm["normalized_qc"] = normalized_qc

    adata.obsm["qc_pc"] = pbmc_metadata["qc_pc"][idx_metadata]
    labels = pbmc_metadata["clusters"][idx_metadata]
    cell_types = pbmc_metadata["list_clusters"]
    adata.obs["labels"] = labels
    adata.uns["cell_types"] = cell_types
    adata.obs["str_labels"] = [cell_types[i] for i in labels]

    adata.var["n_counts"] = np.squeeze(np.asarray(np.sum(adata.X, axis=0)))

    if run_setup_anndata:
        setup_anndata(adata, batch_key="batch", labels_key="labels")
    return adata
