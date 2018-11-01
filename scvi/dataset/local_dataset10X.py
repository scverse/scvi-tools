import os
import pickle
import tarfile

import numpy as np
import pandas as pd
from scipy import io
from scipy.sparse import csr_matrix

from scvi.dataset import GeneExpressionDataset


class LocalDataset10X(GeneExpressionDataset):


    def __init__(self, save_path):
        self.save_path = save_path
        expression_data, gene_names = self.preprocess()
        super(LocalDataset10X, self).__init__(*GeneExpressionDataset.get_attributes_from_matrix(
            expression_data), gene_names=gene_names)


    def preprocess(self):
        print("Preprocessing dataset")
        genes_info = pd.read_csv(self.path + 'genes.tsv', sep='\t', header=None)
        gene_names = genes_info.values[:, 0].astype(np.str).ravel()
        if os.path.exists(self.path + 'barcodes.tsv'):
            self.barcodes = pd.read_csv(self.path + 'barcodes.tsv', sep='\t', header=None)
        self.gene_symbols = genes_info.values[:, 1].astype(np.str).ravel()
        expression_data = io.mmread(self.path + 'matrix.mtx').T
        if self.dense:
            expression_data = expression_data.A
        else:
            expression_data = csr_matrix(expression_data)

        print("Finished preprocessing dataset")
        return expression_data, gene_names
