import collections

import numpy as np
from scipy.sparse import csc_matrix, csr_matrix, vstack
import h5py
from sklearn.preprocessing import StandardScaler
import time
from tqdm import trange

from .dataset import GeneExpressionDataset

batch_idx_10x = [1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0,
                 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0,
                 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0,
                 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1]


def subsample_barcodes(gbm, barcode_indices):
    barcodes = gbm.barcodes[barcode_indices]
    return GeneBCMatrix(gbm.gene_ids, gbm.gene_names, barcodes,
                        gbm.matrix[:, barcode_indices])


def subsample_genes(gbm, genes_indices):
    gene_ids = gbm.gene_ids[genes_indices]
    gene_names = gbm.gene_names[genes_indices]
    return GeneBCMatrix(gene_ids, gene_names, gbm.barcodes,
                        gbm.matrix[genes_indices, :])


def get_expression(gbm, gene_name):
    gene_indices = np.where(gbm.gene_names == gene_name)[0]
    if len(gene_indices) == 0:
        raise Exception("%s was not found in list of gene names." % gene_name)
    return gbm.matrix[gene_indices[0], :].toarray().squeeze()


class BrainLargeDataset(GeneExpressionDataset):
    r""" Loads brain-large dataset.

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

    def __init__(self, subsample_size=None, save_path='data/', nb_genes_kept=720):
        self.subsample_size = 50000#subsample_size
        self.save_path = save_path
        self.nb_genes_kept = nb_genes_kept
        self.url = "http://cf.10xgenomics.com/samples/cell-exp/1.3.0/1M_neurons/" \
                   "1M_neurons_filtered_gene_bc_matrices_h5.h5"
        # originally: "1M_neurons_filtered_gene_bc_matrices_h5.h5"

        self.download_name = "genomics.h5"

        Xs = self.download_and_preprocess()
        super(BrainLargeDataset, self).__init__(
            *GeneExpressionDataset.get_attributes_from_list(Xs)
        )

    def preprocess(self):
        print("Preprocessing Brain Large data")

        filtered_matrix_h5 = self.save_path + self.download_name
        with h5py.File(filtered_matrix_h5) as f:
            dset = f["mm10"]
            n_genes, n_cells = f["mm10"]["shape"]
            if self.subsample_size is None:
                self.subsample_size = n_cells
            indptr = dset['indptr'][...]

            ns_cells = 10000
            ns_indptr = indptr[:(ns_cells + 1)]
            ns_nnz = ns_indptr[-1]
            ns_data = dset["data"][:ns_nnz].astype(np.float32)
            ns_indices = dset["data"][:ns_nnz]
            ns_sparse = csc_matrix((ns_data, ns_indices, ns_indptr), shape=(n_genes, ns_cells))
            ns_dense = ns_sparse.toarray()

            std_scaler = StandardScaler(with_mean=False)
            std_scaler.fit(ns_dense)
            subset_genes = np.argsort(std_scaler.var_)[::-1][:self.nb_genes_kept]

            nb_matrices = []
            nb_cells = 100000
            nb_iters = int(self.subsample_size / nb_cells) + (self.subsample_size % nb_cells > 0)
            for i in range(nb_iters):
                nb_indptr = indptr[(i * nb_cells):((1 + i) * nb_cells + 1)]
                nb_nnz_a = nb_indptr[0]
                nb_nnz_b = nb_indptr[-1]
                nb_indptr = (nb_indptr - nb_nnz_a).astype(np.int32)
                nb2_cells = len(nb_indptr) - 1
                nb_data = dset["data"][nb_nnz_a:nb_nnz_b].astype(np.float32)
                nb_indices = dset["indices"][nb_nnz_a:nb_nnz_b].astype(np.int32)
                nb_sparse = csr_matrix((nb_data, nb_indices, nb_indptr), shape=(nb2_cells, n_genes))
                nb_filtered = nb_sparse[:, subset_genes]
                del nb_sparse
                nb_matrices.append(nb_filtered)
                print("loaded {} / {} cells".format(i * nb_cells + nb2_cells, self.subsample_size))

        matrix = vstack(nb_matrices)

        print("%d cells subsampled" % matrix.shape[0])
        print("%d genes subsampled" % matrix.shape[1])

        return [matrix,]
