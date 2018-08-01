import csv

import numpy as np

from scvi.dataset import GeneExpressionDataset


def reorganize(x, genes, ordered_genes):
    # X must be a numpy matrix
    new_order_first = []
    for ordered_gene in range(len(ordered_genes)):
        for gene in range(len(genes)):
            if ordered_genes[ordered_gene].lower() == genes[gene].lower():
                new_order_first.append(gene)
    new_order_second = [x for x in range(len(genes)) if x not in new_order_first]
    new_order = new_order_first + new_order_second

    return x[:, new_order], genes[new_order]


class CortexDatasetCustom(GeneExpressionDataset):
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

    def __init__(self, save_path='data/', genes_fish=None, genes_to_discard=None,
                 genes_to_keep=None, additional_genes=558):
        # Generating samples according to a ZINB process

        self.save_path = save_path
        self.download_name = 'expression.bin'
        self.url = "https://storage.googleapis.com/linnarsson-lab-www-blobs/blobs/cortex/" \
                   "expression_mRNA_17-Aug-2014.txt"
        self.genes_to_discard = genes_to_discard
        self.genes_to_keep = genes_to_keep
        self.genes_fish = genes_fish
        self.additional_genes = additional_genes
        expression_data, labels, gene_names = self.download_and_preprocess()

        cell_types = ["astrocytes_ependymal",
                      "endothelial-mural",
                      "interneurons",
                      "microglia",
                      "oligodendrocytes",
                      "pyramidalCA1",
                      "pyramidal-SS"]

        super(CortexDatasetCustom, self).__init__(
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

        # Getting the indexes to keep
        selected_fish = []
        if self.genes_fish is not None:
            for gene_fish in self.genes_fish:
                for gene_cortex in range(len(gene_names)):
                    if gene_names[gene_cortex].lower() == gene_fish.lower():
                        selected_fish.append(gene_cortex)
        to_keep = []
        if self.genes_to_keep is not None:
            for gene in self.genes_to_keep:
                for gene_cortex in range(len(gene_names)):
                    if gene_names[gene_cortex].lower() == gene.lower():
                        to_keep.append(gene_cortex)

        selected = np.std(expression_data, axis=0).argsort()[-self.additional_genes:][::-1]
        selected = np.unique(np.concatenate((selected, np.array(to_keep), np.array(selected_fish))))
        selected = np.array([int(select) for select in selected])
        expression_data = expression_data[:, selected]
        gene_names = gene_names[selected]
        expression_data, gene_names = reorganize(expression_data, gene_names, self.genes_fish)
        umi = np.sum(expression_data[:, :len(self.genes_fish)], axis=1)
        expression_data = expression_data[umi > 10, :]
        labels = labels[umi > 10]
        print("Finished preprocessing Cortex data")
        return expression_data, labels, gene_names
