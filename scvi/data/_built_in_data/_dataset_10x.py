import logging
import os
import shutil
import tarfile
import warnings
from typing import Tuple

import numpy as np

from scvi.data._download import _download

logger = logging.getLogger(__name__)

available_datasets = {
    "1.1.0": [
        "frozen_pbmc_donor_a",
        "frozen_pbmc_donor_b",
        "frozen_pbmc_donor_c",
        "fresh_68k_pbmc_donor_a",
        "cd14_monocytes",
        "b_cells",
        "cd34",
        "cd56_nk",
        "cd4_t_helper",
        "regulatory_t",
        "naive_t",
        "memory_t",
        "cytotoxic_t",
        "naive_cytotoxic",
    ],
    "1.3.0": ["1M_neurons"],
    "2.1.0": ["pbmc8k", "pbmc4k", "t_3k", "t_4k", "neuron_9k"],
    "3.0.0": [
        "pbmc_1k_protein_v3",
        "pbmc_10k_protein_v3",
        "malt_10k_protein_v3",
        "pbmc_1k_v2",
        "pbmc_1k_v3",
        "pbmc_10k_v3",
        "hgmm_1k_v2",
        "hgmm_1k_v3",
        "hgmm_5k_v3",
        "hgmm_10k_v3",
        "neuron_1k_v2",
        "neuron_1k_v3",
        "neuron_10k_v3",
        "heart_1k_v2",
        "heart_1k_v3",
        "heart_10k_v3",
    ],
    "3.1.0": ["5k_pbmc_protein_v3", "5k_pbmc_protein_v3_nextgem"],
}

dataset_to_group = {
    dataset_name: group
    for group, list_datasets in available_datasets.items()
    for dataset_name in list_datasets
}

group_to_url_skeleton = {
    "1.1.0": "http://cf.10xgenomics.com/samples/cell-exp/{}/{}/{}_{}_gene_bc_matrices.tar.gz",
    "1.3.0": "http://cf.10xgenomics.com/samples/cell-exp/{}/{}/{}_{}_gene_bc_matrices_h5.h5",
    "2.1.0": "http://cf.10xgenomics.com/samples/cell-exp/{}/{}/{}_{}_gene_bc_matrices.tar.gz",
    "3.0.0": "http://cf.10xgenomics.com/samples/cell-exp/{}/{}/{}_{}_feature_bc_matrix.h5",
    "3.1.0": "http://cf.10xgenomics.com/samples/cell-exp/{}/{}/{}_{}_feature_bc_matrix.h5",
}

group_to_filename_skeleton = {
    "1.1.0": "{}_gene_bc_matrices.tar.gz",
    "1.3.0": "{}_gene_bc_matrices_h5.h5",
    "2.1.0": "{}_gene_bc_matrices.tar.gz",
    "3.0.0": "{}_feature_bc_matrix.h5",
    "3.1.0": "{}_feature_bc_matrix.h5",
}


def _load_dataset_10x(
    dataset_name: str = None,
    filename: str = None,
    save_path: str = "data/10X",
    url: str = None,
    return_filtered: bool = True,
    remove_extracted_data: bool = False,
    **scanpy_read_10x_kwargs,
):
    try:
        import scanpy
    except ImportError:
        raise ImportError("Please install scanpy -- `pip install scanpy`")

    # form data url and filename unless manual override
    if dataset_name is not None:
        if url is not None:
            warnings.warn("dataset_name provided, manual url is disregarded.")
        if filename is not None:
            warnings.warn("dataset_name provided, manual filename is disregarded.")
        group = dataset_to_group[dataset_name]
        url_skeleton = group_to_url_skeleton[group]

        filter_type = "filtered" if return_filtered else "raw"
        url = url_skeleton.format(group, dataset_name, dataset_name, filter_type)
        filename_skeleton = group_to_filename_skeleton[group]
        filename = filename_skeleton.format(filter_type)
        save_path = os.path.join(save_path, dataset_name)
    elif filename is not None and url is not None:
        logger.info("Loading 10X dataset with custom url and filename")
    elif filename is not None and url is None:
        logger.info("Loading local 10X dataset with custom filename")
    else:
        logger.info("Loading extracted local 10X dataset with custom filename")
    _download(url, save_path=save_path, filename=filename)
    file_path = os.path.join(save_path, filename)

    # untar
    download_is_targz = url[-7:] == ".tar.gz"
    was_extracted = False
    if download_is_targz is True:
        if not os.path.exists(file_path[:-7]):  # nothing extracted yet
            if tarfile.is_tarfile(file_path):
                logger.info("Extracting tar file")
                tar = tarfile.open(file_path, "r:gz")
                tar.extractall(path=save_path)
                was_extracted = True
                tar.close()
        path_to_data_folder, suffix = _find_path_to_mtx(save_path)
        adata = scanpy.read_10x_mtx(path_to_data_folder, **scanpy_read_10x_kwargs)
        if was_extracted and remove_extracted_data:
            folders_in_save_path = path_to_data_folder[len(save_path) + 1 :].split("/")
            extracted_folder_path = save_path + "/" + folders_in_save_path[0]
            logger.info(f"Removing extracted data at {extracted_folder_path}")
            shutil.rmtree(extracted_folder_path)
    else:
        adata = scanpy.read_10x_h5(file_path, **scanpy_read_10x_kwargs)

    adata.var_names_make_unique()
    scanpy.pp.filter_cells(adata, min_counts=1)
    scanpy.pp.filter_genes(adata, min_counts=1)

    return adata


def _find_path_to_mtx(save_path: str) -> Tuple[str, str]:
    """
    Returns exact path for the data in the archive.

    This is required because 10X doesn't have a consistent way of storing their data.
    Additionally, the function returns whether the data is stored in compressed format.

    Returns
    -------
    path in which files are contains and their suffix if compressed.

    """
    for root, _, files in os.walk(save_path):
        # do not consider hidden files
        files = [f for f in files if not f[0] == "."]
        contains_mat = [
            filename == "matrix.mtx" or filename == "matrix.mtx.gz"
            for filename in files
        ]
        contains_mat = np.asarray(contains_mat).any()
        if contains_mat:
            is_tar = files[0][-3:] == ".gz"
            suffix = ".gz" if is_tar else ""
            return root, suffix
    raise FileNotFoundError("No matrix.mtx(.gz) found in path (%s)." % save_path)
