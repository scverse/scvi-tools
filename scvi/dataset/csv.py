from .dataset import GeneExpressionDataset
import pandas as pd
import numpy as np


class CsvDataset(GeneExpressionDataset):

    def __init__(self, filename, save_path='data/', url=None, new_n_genes=600, subset_genes=None, compression=None):
        self.download_name = filename  # The given csv file is
        self.save_path = save_path
        self.url = url
        self.compression = compression

        data, gene_names = self.download_and_preprocess()

        super(CsvDataset, self).__init__(
            *GeneExpressionDataset.get_attributes_from_matrix(
                data), gene_names=gene_names)

        self.subsample_genes(new_n_genes, subset_genes)

    def preprocess(self):
        print("Preprocessing dataset")

        data = pd.read_csv(self.save_path + self.download_name, index_col=0, compression=self.compression).T
        gene_names = np.array(data.columns, dtype=str)

        data = data.values
        print("Finished preprocessing dataset")
        return data, gene_names
