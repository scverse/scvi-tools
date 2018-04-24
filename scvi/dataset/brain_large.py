import collections
import os
import time

import numpy as np
import scipy.sparse as sp_sparse
import tables

from .dataset import GeneExpressionDataset

GeneBCMatrix = collections.namedtuple('GeneBCMatrix', ['gene_ids', 'gene_names', 'barcodes', 'matrix'])


def get_matrix_from_h5(filename, genome):
    with tables.open_file(filename, 'r') as f:
        try:
            dsets = {}
            for node in f.walk_nodes('/' + genome, 'Array'):
                dsets[node.name] = node.read()
            matrix = sp_sparse.csc_matrix((dsets['data'], dsets['indices'], dsets['indptr']),
                                          shape=dsets['shape'])
            return GeneBCMatrix(dsets['genes'], dsets['gene_names'], dsets['barcodes'], matrix)
        except tables.NoSuchNodeError:
            raise Exception("Genome %s does not exist in this file." % genome)
        except KeyError:
            raise Exception("File is missing one or more required datasets.")


def subsample_barcodes(gbm, barcode_indices, test=False):
    barcodes = gbm.barcodes[barcode_indices] if not test else gbm.barcodes
    return GeneBCMatrix(gbm.gene_ids, gbm.gene_names, barcodes,
                        gbm.matrix[:, barcode_indices])


def get_expression(gbm, gene_name):
    gene_indices = np.where(gbm.gene_names == gene_name)[0]
    if len(gene_indices) == 0:
        raise Exception("%s was not found in list of gene names." % gene_name)
    return gbm.matrix[gene_indices[0], :].toarray().squeeze()


def save_matrix_to_h5(gbm, filename, genome):
    flt = tables.Filters(complevel=1)
    with tables.open_file(filename, 'w', filters=flt) as f:
        try:
            group = f.create_group(f.root, genome)
            f.create_carray(group, 'genes', obj=gbm.gene_ids)
            f.create_carray(group, 'gene_names', obj=gbm.gene_names)
            f.create_carray(group, 'barcodes', obj=gbm.barcodes)
            f.create_carray(group, 'data', obj=gbm.matrix.data)
            f.create_carray(group, 'indices', obj=gbm.matrix.indices)
            f.create_carray(group, 'indptr', obj=gbm.matrix.indptr)
            f.create_carray(group, 'shape', obj=gbm.matrix.shape)
        except ValueError:
            raise Exception("Failed to write H5 file.")


class BrainLargeDataset(GeneExpressionDataset):
    url = "http://cf.10xgenomics.com/samples/cell-exp/1.3.0/1M_neurons/1M_neurons_filtered_gene_bc_matrices_h5.h5"

    def __init__(self, subsample_size=None, test=False):
        """
        :param subsample_size: In thousands of barcodes kept (by default 1*1000=1000 kept)
        :param test: A boolean to indicate if we use pytest subsampled file
        """
        self.subsample_size = subsample_size if not test else 128
        self.save_path = 'data/'
        self.test = test
        # originally: "1M_neurons_filtered_gene_bc_matrices_h5.h5"

        if not self.test:
            self.download_name = "genomics.h5"
        else:
            self.download_name = "../tests/data/genomics_subsampled.h5"

        if subsample_size is None:
            self.final_name = self.download_name
        else:
            self.final_name = "genomics_subsampled_%d.h5" % subsample_size

        self.genome = "mm10"
        self.download_and_preprocess()
        h5_object = get_matrix_from_h5(self.save_path + self.final_name, self.genome)
        super(BrainLargeDataset, self).__init__(
            *GeneExpressionDataset.get_attributes_from_matrix(h5_object.matrix.transpose().tocsr(copy=False))
        )

    def preprocess(self):
        print("Preprocessing Brain Large data")
        tic = time.time()
        np.random.seed(0)

        filtered_matrix_h5 = self.save_path + self.download_name

        gene_bc_matrix = get_matrix_from_h5(filtered_matrix_h5, self.genome)

        subsample_bcs = self.subsample_size
        subset = np.sort(np.random.choice(gene_bc_matrix.matrix.shape[1], size=subsample_bcs, replace=self.test))
        subsampled_matrix = subsample_barcodes(gene_bc_matrix, subset, test=self.test)

        save_matrix_to_h5(subsampled_matrix, self.save_path + self.final_name, "mm10")
        toc = time.time()
        print("Preprocessing finished in : %d sec." % int(toc - tic))

    def download_and_preprocess(self):
        if not os.path.exists(self.save_path + self.final_name):
            if not os.path.exists(self.save_path + self.download_name):
                self.download()
            self.preprocess()
