import logging
import os

import anndata
import h5py
import numpy as np
import scipy.sparse as sp_sparse

from scvi.data._download import _download

logger = logging.getLogger(__name__)


def _load_brainlarge_dataset(
    save_path: str = "data/",
    sample_size_gene_var: int = 10000,
    max_cells_to_keep: int = None,
    n_genes_to_keep: int = 720,
    loading_batch_size: int = 100000,
) -> anndata.AnnData:
    """Loads brain-large dataset."""
    url = "http://cf.10xgenomics.com/samples/cell-exp/1.3.0/1M_neurons/1M_neurons_filtered_gene_bc_matrices_h5.h5"
    save_fn = "brain_large.h5"

    _download(url, save_path, save_fn)
    adata = _load_brainlarge_file(
        os.path.join(save_path, save_fn),
        sample_size_gene_var=sample_size_gene_var,
        max_cells_to_keep=max_cells_to_keep,
        n_genes_to_keep=n_genes_to_keep,
        loading_batch_size=loading_batch_size,
    )
    return adata


def _load_brainlarge_file(
    path_to_file: str,
    sample_size_gene_var: int,
    max_cells_to_keep: int,
    n_genes_to_keep: int,
    loading_batch_size: int,
) -> anndata.AnnData:
    logger.info("Preprocessing Brain Large data")
    print(path_to_file)
    with h5py.File(path_to_file, "r") as f:
        data = f["mm10"]
        nb_genes, nb_cells = f["mm10"]["shape"]
        n_cells_to_keep = (
            max_cells_to_keep if max_cells_to_keep is not None else nb_cells
        )
        index_partitioner = data["indptr"][...]
        # estimate gene variance using a subset of cells.
        index_partitioner_gene_var = index_partitioner[: (sample_size_gene_var + 1)]
        last_index_gene_var_sample = index_partitioner_gene_var[-1]
        gene_var_sample_matrix = sp_sparse.csc_matrix(
            (
                data["data"][:last_index_gene_var_sample].astype(np.float32),
                data["indices"][:last_index_gene_var_sample],
                index_partitioner_gene_var,
            ),
            shape=(nb_genes, len(index_partitioner_gene_var) - 1),
        )
        mean = gene_var_sample_matrix.mean(axis=1)
        var = gene_var_sample_matrix.multiply(gene_var_sample_matrix).mean(
            axis=1
        ) - np.multiply(mean, mean)
        subset_genes = np.squeeze(np.asarray(var)).argsort()[-n_genes_to_keep:][::-1]
        del gene_var_sample_matrix, mean, var

        n_iters = int(n_cells_to_keep / loading_batch_size) + (
            n_cells_to_keep % loading_batch_size > 0
        )
        for i in range(n_iters):
            index_partitioner_batch = index_partitioner[
                (i * loading_batch_size) : ((1 + i) * loading_batch_size + 1)
            ]
            first_index_batch = index_partitioner_batch[0]
            last_index_batch = index_partitioner_batch[-1]
            index_partitioner_batch = (
                index_partitioner_batch - first_index_batch
            ).astype(np.int32)
            n_cells_batch = len(index_partitioner_batch) - 1
            data_batch = data["data"][first_index_batch:last_index_batch].astype(
                np.float32
            )
            indices_batch = data["indices"][first_index_batch:last_index_batch].astype(
                np.int32
            )
            matrix_batch = sp_sparse.csr_matrix(
                (data_batch, indices_batch, index_partitioner_batch),
                shape=(n_cells_batch, nb_genes),
            )[:, subset_genes]
            # stack on the fly to limit RAM usage
            if i == 0:
                matrix = matrix_batch
            else:
                matrix = sp_sparse.vstack([matrix, matrix_batch])
            logger.info(
                "loaded {} / {} cells".format(
                    i * loading_batch_size + n_cells_batch, n_cells_to_keep
                )
            )
    logger.info("%d cells subsampled" % matrix.shape[0])
    logger.info("%d genes subsampled" % matrix.shape[1])
    adata = anndata.AnnData(matrix)
    adata.obs["labels"] = np.zeros(matrix.shape[0])
    adata.obs["batch"] = np.zeros(matrix.shape[0])

    counts = adata.X.sum(1)
    if sp_sparse.issparse(counts):
        counts = counts.A1

    gene_num = (adata.X > 0).sum(1)
    if sp_sparse.issparse(gene_num):
        gene_num = gene_num.A1

    adata = adata[counts > 1]
    adata = adata[gene_num > 1]

    return adata.copy()
