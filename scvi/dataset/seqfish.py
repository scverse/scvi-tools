import pandas as pd
import os
import logging

from scvi.dataset.dataset import GeneExpressionDataset


class SeqfishDataset(GeneExpressionDataset):

    def __init__(self, save_path='data/'):
        self.download_name = "SeqFISH.xlsx"
        self.save_path = save_path
        self.url = "https://www.cell.com/cms/attachment/2080562255/2072099886/mmc6.xlsx"

        data = self.download_and_preprocess()

        super().__init__(*GeneExpressionDataset.get_attributes_from_matrix(data))

    def preprocess(self):
        logging.info("Preprocessing dataset")

        xl = pd.ExcelFile(os.path.join(self.save_path, self.download_name))
        ds = xl.parse("Hippocampus Counts")  # They also used cell by genes

        logging.info("Finished preprocessing dataset")
        return ds.values[:, 1:].astype(int)
