# On the 10X website:
# The main categories (Cell Ranger 1.1.0 / Cell Ranger 2.1.0 / ...)
# have same access suffix for each of their dataset.
# For dataset name (eg. 'pbmc8k', 'pbmc4k', ect...) their are two available specifications,
# either filtered or raw data
import os
import tarfile

import numpy as np
import pandas as pd
from scipy.io import mmread
from scipy.sparse import csr_matrix

from scvi.dataset import GeneExpressionDataset

available_datasets = {"1.1.0":
                      ["frozen_pbmc_donor_a",
                       "frozen_pbmc_donor_b",
                       "frozen_pbmc_donor_c"],
                      "2.1.0":
                      ["pbmc8k",
                       "pbmc4k",
                       "t_3k",
                       "t_4k",
                       "neuron_9k"]}

to_groups = dict([(dataset_name, group)
                  for group, list_datasets in available_datasets.items()
                  for dataset_name in list_datasets])
available_specification = ['filtered', 'raw']


class Dataset10X(GeneExpressionDataset):
    def __init__(self, name, type='filtered'):
        group = to_groups[name]
        self.url = ("http://cf.10xgenomics.com/samples/cell-exp/%s/%s/%s_%s_gene_bc_matrices.tar.gz" %
                    (group, name, name, type))
        self.save_path = 'data/10X/%s/' % name
        self.download_name = 'gene_bc_matrices.tar.gz'
        self.save_name = 'gene_bc_matrices'
        expression_data, gene_names = self.download_and_preprocess()
        super(Dataset10X, self).__init__(*GeneExpressionDataset.get_attributes_from_matrix(
            expression_data), gene_names=gene_names)

    def preprocess(self):
        path = self.save_path + self.download_name
        if len(os.listdir(self.save_path)) == 1:  # nothing extracted yet
            print("Extracting tar file")
            tar = tarfile.open(path, "r:gz")
            tar.extractall()
            tar.close()
        path = self.save_path + [dir for dir in os.listdir(self.save_path) if dir != self.download_name][0] + '/'
        path += os.listdir(path)[0] + '/'
        genes_info = pd.read_csv(path + 'genes.tsv', sep='\t', header=None)
        gene_names = genes_info.values[:, 0].astype(np.str).ravel()
        expression_data = csr_matrix(mmread(path + "matrix.mtx").T)
        return expression_data, gene_names
