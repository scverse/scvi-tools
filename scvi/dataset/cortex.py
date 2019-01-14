import csv
import os
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

    def __init__(self, save_path='data/', genes_to_keep=[], genes_fish=[], additional_genes=558):
        # Generating samples according to a ZINB process
        self.save_path = save_path
        self.download_name = 'expression.bin'
        self.url = "https://storage.googleapis.com/linnarsson-lab-www-blobs/blobs/cortex/" \
                   "expression_mRNA_17-Aug-2014.txt"
        # If we want to harmonize the dataset with the OsmFISH dataset, we need to keep
        # OsmFISH genes and order the genes from Cortex accordingly
        self.genes_fish = genes_fish
        # If there are specific genes we'd like to keep
        self.genes_to_keep = genes_to_keep
        # Number of genes we want to keep
        self.additional_genes = additional_genes
        expression_data, labels, gene_names, cell_types = self.download_and_preprocess()

        super().__init__(
            *GeneExpressionDataset.get_attributes_from_matrix(
                expression_data,
                labels=labels),
            gene_names=np.char.upper(gene_names), cell_types=cell_types)

    def preprocess(self):
        print("Preprocessing Cortex data")
        rows = []
        gene_names = []
        with open(os.path.join(self.save_path, self.download_name), 'r') as csvfile:
            data_reader = csv.reader(csvfile, delimiter='\t')
            clusters = None
            for i, row in enumerate(data_reader):
                if i == 1:
                    precise_clusters = np.array(row, dtype=str)[2:]
                if i == 8:  # 7 + 1 in pandas
                    clusters = np.array(row, dtype=str)[2:]
                if i >= 11:  # 10 + 1 in pandas
                    rows.append(row[1:])
                    gene_names.append(row[0])
        cell_types, labels = np.unique(clusters, return_inverse=True)
        _, self.precise_labels = np.unique(precise_clusters, return_inverse=True)

        expression_data = np.array(rows, dtype=np.int).T[1:]
        gene_names = np.array(gene_names, dtype=np.str)

        additional_genes = []
        for gene_cortex in range(len(gene_names)):
            for gene_fish in self.genes_fish:
                if gene_names[gene_cortex].lower() == gene_fish.lower():
                    additional_genes.append(gene_cortex)
        for gene_cortex in range(len(gene_names)):
            for gene_fish in self.genes_to_keep:
                if gene_names[gene_cortex].lower() == gene_fish.lower():
                    additional_genes.append(gene_cortex)

        selected = np.std(expression_data, axis=0).argsort()[-self.additional_genes:][::-1]
        selected = np.unique(np.concatenate((selected, np.array(additional_genes))))
        selected = np.array([int(select) for select in selected])
        expression_data = expression_data[:, selected]
        gene_names = gene_names[selected]

        # Then we reorganize the genes so that the genes from the smFISH dataset
        # appear first
        if len(self.genes_fish) > 0:
            expression_data, gene_names = self.reorder_genes(expression_data, gene_names, self.genes_fish)
            umi = np.sum(expression_data[:, :len(self.genes_fish)], axis=1)
            expression_data = expression_data[umi > 10, :]
            labels = labels[umi > 10]

        print("Finished preprocessing Cortex data")
        return expression_data, labels, gene_names, cell_types

    @staticmethod
    def reorder_genes(x, genes, first_genes):
        """
        In case the order of the genes needs to be changed:
        puts the gene present in ordered_genes first, conserving
        the same order.
        """
        # X must be a numpy matrix
        new_order_first = []
        for ordered_gene in range(len(first_genes)):
            for gene in range(len(genes)):
                if first_genes[ordered_gene].lower() == genes[gene].lower():
                    new_order_first.append(gene)
        new_order_second = [x for x in range(len(genes)) if x not in new_order_first]
        new_order = new_order_first + new_order_second

        return x[:, new_order], genes[new_order]
