import numpy as np

import tarfile
from scipy import io

from .dataset import GeneExpressionDataset


class BrainSmallDataset(GeneExpressionDataset):
    url = "http://cf.10xgenomics.com/samples/cell-exp/2.1.0/neuron_9k/" + \
          "neuron_9k_filtered_gene_bc_matrices.tar.gz"

    def __init__(self, save_path='data/'):
        self.save_path = save_path

        fn = "filtered_gene_bc_matrices"
        self.download_name = fn + ".tar.gz"
        self.gene_filename = fn + "/mm10/genes.tsv"
        self.expression_filename = fn + "/mm10/matrix.mtx"

        expression_data, gene_names = self.download_and_preprocess()

        super(BrainSmallDataset, self).__init__(
            *GeneExpressionDataset.get_attributes_from_matrix(
                expression_data), gene_names=gene_names)

    def preprocess(self):
        print("Preprocessing Brain Small data")
        gene_names = []

        tar = tarfile.open(self.save_path + self.download_name)
        gene_file = tar.extractfile(tar.getmember(self.gene_filename))

        lines = gene_file.readlines()
        for line in lines:
            gene_names.append(line[:-2])  # remove "\n"

        expression_file = tar.extractfile(tar.getmember(self.expression_filename))
        expression_data = io.mmread(expression_file).T.A

        gene_names = np.array(gene_names, dtype=np.str)

        selected = np.std(expression_data, axis=0).argsort()[-3000:][::-1]
        expression_data = expression_data[:, selected]
        gene_names = gene_names[selected]

        return expression_data, gene_names
