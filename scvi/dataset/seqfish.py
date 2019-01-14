import pandas as pd
from .dataset import GeneExpressionDataset
import os


class SeqfishDataset(GeneExpressionDataset):

    def __init__(self, save_path='data/'):
        self.download_name = "SeqFISH.xlsx"
        self.save_path = save_path
        self.url = "https://www.cell.com/cms/attachment/2080562255/2072099886/mmc6.xlsx"

        data = self.download_and_preprocess()

        super().__init__(*GeneExpressionDataset.get_attributes_from_matrix(data))

    def preprocess(self):
        print("Preprocessing dataset")

        xl = pd.ExcelFile(os.path.join(self.save_path, self.download_name))
        ds = xl.parse("Hippocampus Counts")  # They also used cell by genes

        print("Finished preprocessing dataset")
        return ds.values[:, 1:].astype(int)
