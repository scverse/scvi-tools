import pandas as pd
import numpy as np

from .dataset import GeneExpressionDataset


class CbmcDataset(GeneExpressionDataset):
    url = "https://www.ncbi.nlm.nih.gov/geo/download/" + \
          "?acc=GSE100866&format=file&file=GSE100866%5FCBMC%5F8K%5F13AB%5F10X%2DRNA%5Fumi%2Ecsv%2Egz"

    def __init__(self, save_path='data/'):
        self.save_path = save_path

        # originally: GSE100866_CBMC_8K_13AB_10X-RNA_umi.csv.gz

        self.download_name = 'cbmc.csv.gz'

        expression_data, gene_names = self.download_and_preprocess()

        super(CbmcDataset, self).__init__(
            *GeneExpressionDataset.get_attributes_from_matrix(
                expression_data), gene_names=gene_names)

    def preprocess(self):
        print("Preprocessing cbmc data")
        expression = pd.read_csv(self.save_path + self.download_name, index_col=0, compression='gzip').T
        gene_names = np.array(expression.columns, dtype=str)

        selected = np.std(expression.values, axis=0).argsort()[-600:][::-1]
        expression_data = expression.values[:, selected]

        return expression_data, gene_names
