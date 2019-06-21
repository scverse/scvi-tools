import logging
import os

import h5py
import numpy as np
import scipy.sparse as sp_sparse

from scvi.dataset.dataset import DownloadableDataset

logger = logging.getLogger(__name__)


class BrainLargeDataset(DownloadableDataset):
    """Loads brain-large dataset.

    This dataset contains 1.3 million brain cells from `10x Genomics`_. We randomly shuffle the data to get a 1M
    subset of cells and order genes by variance to retain first 10,000 and then 720 sampled variable genes. This
    dataset is then sampled multiple times in cells for the runtime and goodness-of-fit analysis. We report imputation
    scores on the 10k cells and 720 genes samples only.

    Args:
        :save_path: Save path of raw data file. Default: ``'data/'``.

    Examples:
        >>> gene_dataset = BrainLargeDataset()

    .. _10x Genomics:
        https://support.10xgenomics.com/single-cell-gene-expression/datasets

    """

    def __init__(
        self,
        save_path: str = "data/",
        filename: str = None,
        delayed_populating: bool = False,
        sample_size_gene_var: int = 10000,
        max_cells_to_keep: int = None,
        nb_genes_to_keep: int = 720,
        loading_batch_size: int = 100000,
    ):
        # used in populate, should not be moved after the call to super().__init__()
        self.sample_size_gene_var = sample_size_gene_var
        self.max_cells_to_keep = max_cells_to_keep
        self.nb_genes_to_keep = nb_genes_to_keep
        self.loading_batch_size = loading_batch_size
        super().__init__(
            urls=(
                "http://cf.10xgenomics.com/samples/cell-exp/1.3.0/1M_neurons/"
                "1M_neurons_filtered_gene_bc_matrices_h5.h5"
            ),
            filenames=filename if filename is not None else "brain_large.h5",
            save_path=save_path,
            delayed_populating=delayed_populating,
        )

    def populate(self):
        logger.info("Preprocessing Brain Large data")
        with h5py.File(os.path.join(self.save_path, self.filenames[0]), "r") as f:
            data = f["mm10"]
            nb_genes, nb_cells = f["mm10"]["shape"]
            self.n_cells_to_keep = (
                self.max_cells_to_keep
                if self.max_cells_to_keep is not None
                else nb_cells
            )
            index_partitioner = data["indptr"][...]
            # estimate gene variance using a subset of cells.
            index_partitioner_gene_var = index_partitioner[
                : (self.sample_size_gene_var + 1)
            ]
            last_index_gene_var_sample = index_partitioner_gene_var[-1]
            gene_var_sample_matrix = sp_sparse.csc_matrix(
                (
                    data["data"][:last_index_gene_var_sample].astype(np.float32),
                    data["indices"][:last_index_gene_var_sample],
                    index_partitioner_gene_var,
                ),
                shape=(nb_genes, self.sample_size_gene_var),
            )
            mean = gene_var_sample_matrix.mean(axis=1)
            var = gene_var_sample_matrix.multiply(gene_var_sample_matrix).mean(
                axis=1
            ) - np.multiply(mean, mean)
            self.subset_genes = (
                np.squeeze(np.asarray(var)).argsort()[-self.nb_genes_to_keep :][::-1]
            )
            del gene_var_sample_matrix, mean, var

            n_iters = int(self.n_cells_to_keep / self.loading_batch_size) + (
                self.n_cells_to_keep % self.loading_batch_size > 0
            )
            for i in range(n_iters):
                index_partitioner_batch = index_partitioner[
                    (i * self.loading_batch_size) : (
                        (1 + i) * self.loading_batch_size + 1
                    )
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
                indices_batch = data["indices"][
                    first_index_batch:last_index_batch
                ].astype(np.int32)
                matrix_batch = sp_sparse.csr_matrix(
                    (data_batch, indices_batch, index_partitioner_batch),
                    shape=(n_cells_batch, nb_genes),
                )[:, self.subset_genes]
                # stack on the fly to limit RAM usage
                if i == 0:
                    matrix = matrix_batch
                else:
                    matrix = sp_sparse.vstack([matrix, matrix_batch])
                logger.info(
                    "loaded {} / {} cells".format(
                        i * self.loading_batch_size + n_cells_batch,
                        self.n_cells_to_keep,
                    )
                )
        logger.info("%d cells subsampled" % matrix.shape[0])
        logger.info("%d genes subsampled" % matrix.shape[1])
        self.populate_from_data(matrix)
        self.filter_cells_by_count()
