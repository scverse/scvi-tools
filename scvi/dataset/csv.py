from .dataset import GeneExpressionDataset
import pandas as pd
import numpy as np


class CsvDataset(GeneExpressionDataset):

    def __init__(self, filename, save_path='data/', url=None,
                 new_n_genes=600, subset_genes=None, sep=',', compression=None, gene_by_cell=True):
        self.download_name = filename  # The given csv file is
        self.save_path = save_path
        self.url = url
        self.compression = compression
        self.sep = sep
        self.gene_by_cell = gene_by_cell  # Whether the original dataset is genes by cells

        data, gene_names = self.download_and_preprocess()

        super(CsvDataset, self).__init__(
            *GeneExpressionDataset.get_attributes_from_matrix(
                data), gene_names=gene_names)

        self.subsample_genes(new_n_genes, subset_genes)

    def preprocess(self):
        print("Preprocessing dataset")

        if self.gene_by_cell:
            data = pd.read_csv(self.save_path + self.download_name,
                               sep=self.sep, index_col=0, compression=self.compression).T
        else:
            data = pd.read_csv(self.save_path + self.download_name,
                               sep=self.sep, index_col=0, compression=self.compression)

        gene_names = np.array(data.columns, dtype=str)

        data = data.values
        print("Finished preprocessing dataset")
        return data, gene_names


class BreastCancerDataset(CsvDataset):
    def __init__(self, save_path='data/'):
        super(BreastCancerDataset, self).__init__("Layer2_BC_count_matrix-1.tsv", save_path=save_path,
                                                  url="http://www.spatialtranscriptomicsresearch.org/wp-content/"
                                                      "uploads/2016/07/Layer2_BC_count_matrix-1.tsv",
                                                  sep='\t', gene_by_cell=False)


class MouseOBDataset(CsvDataset):
    def __init__(self, save_path='data/'):
        super(MouseOBDataset, self).__init__("Rep11_MOB_count_matrix-1.tsv", save_path=save_path,
                                             url="http://www.spatialtranscriptomicsresearch.org/wp-content/uploads/"
                                                 "2016/07/Rep11_MOB_count_matrix-1.tsv",
                                             sep='\t', gene_by_cell=False)
