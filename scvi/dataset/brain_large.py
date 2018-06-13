import collections
import time

import numpy as np
from scipy.sparse import csc_matrix
import h5py
from sklearn.preprocessing import StandardScaler

from .dataset import GeneExpressionDataset

GeneBCMatrix = collections.namedtuple('GeneBCMatrix', ['gene_ids', 'gene_names', 'barcodes', 'matrix'])

batch_idx_10x = [1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0,
                 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0,
                 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0,
                 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1]


def subsample_barcodes(gbm, barcode_indices, unit_test=False):
    barcodes = gbm.barcodes[barcode_indices] if not unit_test else gbm.barcodes
    return GeneBCMatrix(gbm.gene_ids, gbm.gene_names, barcodes,
                        gbm.matrix[:, barcode_indices])


def subsample_genes(gbm, genes_indices, unit_test=False):
    gene_ids = gbm.gene_ids[genes_indices] if not unit_test else gbm.gene_ids
    gene_names = gbm.gene_names[genes_indices] if not unit_test else gbm.gene_names
    return GeneBCMatrix(gene_ids, gene_names, gbm.barcodes,
                        gbm.matrix[genes_indices, :])


def get_expression(gbm, gene_name):
    gene_indices = np.where(gbm.gene_names == gene_name)[0]
    if len(gene_indices) == 0:
        raise Exception("%s was not found in list of gene names." % gene_name)
    return gbm.matrix[gene_indices[0], :].toarray().squeeze()


class BrainLargeDataset(GeneExpressionDataset):
    url = "http://cf.10xgenomics.com/samples/cell-exp/1.3.0/1M_neurons/1M_neurons_filtered_gene_bc_matrices_h5.h5"

    def __init__(self, subsample_size=None, save_path='data/', nb_genes_kept=720):
        """
        :param subsample_size: In thousands of barcodes kept (by default 1*1000=1000 kept)
        :param unit_test: A boolean to indicate if we use pytest subsampled file
        """
        self.subsample_size = subsample_size
        self.save_path = save_path
        self.nb_genes_kept = nb_genes_kept
        # originally: "1M_neurons_filtered_gene_bc_matrices_h5.h5"

        self.download_name = "genomics.h5"

        Xs = self.download_and_preprocess()
        super(BrainLargeDataset, self).__init__(
            *GeneExpressionDataset.get_attributes_from_list(Xs)
        )

    def preprocess(self):
        print("Preprocessing Brain Large data")
        tic = time.time()

        filtered_matrix_h5 = self.save_path + self.download_name
        with h5py.File(filtered_matrix_h5) as f:
            dset = f["mm10"]
            matrix = csc_matrix((dset['data'], dset['indices'], dset['indptr']), shape=dset['shape'])
            gene_bc_matrix = GeneBCMatrix(dset['genes'][:], dset['gene_names'][:], dset['barcodes'][:], matrix)

        # Downsample *barcodes* from 1306127 to 100000 (~1/10) to
        matrix = gene_bc_matrix.matrix[:, :100000]
        std_scaler = StandardScaler(with_mean=False)
        std_scaler.fit(matrix.transpose().astype(np.float64))
        subset_genes = np.argsort(std_scaler.var_)[::-1][:self.nb_genes_kept]
        subsampled_matrix = subsample_genes(gene_bc_matrix, subset_genes, unit_test=(self.save_path == 'tests/data/'))
        print("%d genes subsampled" % subsampled_matrix.matrix.shape[0])
        print("%d cells subsampled" % subsampled_matrix.matrix.shape[1])

        batch_id = np.random.randint(0, 2, size=subsampled_matrix.matrix.T.shape[0])

        # Choose if Xs should be sparse
        if self.subsample_size is not None:
            idx = np.random.permutation(len(batch_id))[:self.subsample_size]
            X = subsampled_matrix.matrix.T[idx]
            b = batch_id[idx]
            Xs = [X[b == batch].A for batch in (0, 1)]
        else:
            X = subsampled_matrix.matrix.T
            b = batch_id
            Xs = [X[b == batch] for batch in (0, 1)]

        toc = time.time()
        print("Preprocessing finished in : %d sec." % int(toc - tic))

        return Xs  # the list of X (n_barcodes*n_genes)
