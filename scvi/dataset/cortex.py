import csv

import numpy as np

from .dataset import GeneExpressionDataset


class CortexDataset(GeneExpressionDataset):
    r""" Loads cortex dataset.

    The `Mouse Cortex Cells dataset`_ contains 3005 mouse cortex cells and gold-standard labels for
    seven distinct cell types. Each cell type corresponds to a cluster to recover. We retain top 558 genes
    ordered by variance.

    Args:
        :save_path: Save path of raw data file. Default: ``'data/'``.

    Examples:
        >>> gene_dataset = CortexDataset()

    .. _Mouse Cortex Cells dataset:
        https://storage.googleapis.com/linnarsson-lab-www-blobs/blobs/cortex/expression_mRNA_17-Aug-2014.txt

    """

    def __init__(self, save_path='data/'):
        # Generating samples according to a ZINB process
        self.save_path = save_path
        self.download_name = 'expression.bin'
        self.url = "https://storage.googleapis.com/linnarsson-lab-www-blobs/blobs/cortex/" \
                   "expression_mRNA_17-Aug-2014.txt"
        expression_data, labels, gene_names = self.download_and_preprocess()

        cell_types = ["astrocytes_ependymal",
                      "endothelial-mural",
                      "interneurons",
                      "microglia",
                      "oligodendrocytes",
                      "pyramidalCA1",
                      "pyramidal-SS"]

        super(CortexDataset, self).__init__(
            *GeneExpressionDataset.get_attributes_from_matrix(
                expression_data,
                labels=labels),
            gene_names=np.char.upper(gene_names), cell_types=cell_types)

    def preprocess(self):
        print("Preprocessing Cortex data")
        rows = []
        gene_names = []
        with open(self.save_path + self.download_name, 'r') as csvfile:
            data_reader = csv.reader(csvfile, delimiter='\t')
            clusters = None
            for i, row in enumerate(data_reader):
                if i == 8:  # 7 + 1 in pandas
                    clusters = np.array(row, dtype=str)[2:]
                if i >= 11:  # 10 + 1 in pandas
                    rows.append(row[1:])
                    gene_names.append(row[0])

        cell_types, labels = np.unique(clusters, return_inverse=True)

        expression_data = np.array(rows, dtype=np.int).T[1:]
        gene_names = np.array(gene_names, dtype=np.str)

        selected = np.std(expression_data, axis=0).argsort()[-558:][::-1]
        expression_data = expression_data[:, selected]
        gene_names = gene_names[selected]

        print("Finished preprocessing Cortex data")
        return expression_data, labels, gene_names
