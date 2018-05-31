import csv
import gzip

import numpy as np

from .dataset import GeneExpressionDataset


class CbmcDataset(GeneExpressionDataset):
    url = "https://www.ncbi.nlm.nih.gov/geo/download/" + \
          "?acc=GSE100866&format=file&file=GSE100866%5FCBMC%5F8K%5F13AB%5F10X%2DRNA%5Fumi%2Ecsv%2Egz"

    def __init__(self, unit_test=False):
        self.save_path = 'data/'
        self.unit_test = unit_test

        # originally: GSE100866_CBMC_8K_13AB_10X-RNA_umi.csv.gz

        if not self.unit_test:
            self.download_name = 'cbmc.csv.gz'
        else:
            self.download_name = "../tests/data/cbmc_subsampled.csv.gz"

        expression_data, gene_names = self.download_and_preprocess()

        super(CbmcDataset, self).__init__(
            *GeneExpressionDataset.get_attributes_from_matrix(
                expression_data), gene_names=gene_names)

    def preprocess(self):
        print("Preprocessing cbmc data")
        rows = []
        gene_names = []
        with gzip.open(self.save_path + self.download_name, "rt", encoding="utf8") as csvfile:
            data_reader = csv.reader(csvfile, delimiter=',')
            for row in data_reader:
                rows.append(row[1:])
                gene_names.append(row[0])

        expression_data = np.array(rows[1:], dtype=np.int).T
        gene_names = np.array(gene_names[1:], dtype=np.str)

        selected = np.std(expression_data, axis=0).argsort()[-600:][::-1]
        expression_data = expression_data[:, selected]
        gene_names = gene_names[selected]

        return expression_data, gene_names
