from .dataset import GeneExpressionDataset
import anndata
import os
import zipfile
import numpy as np


class AnnDataset(GeneExpressionDataset):

    def __init__(self, filename, download_name=None, save_path='data/', url=None):
        self.download_name = download_name
        self.filename = filename
        self.save_path = save_path
        self.url = url

        data, gene_names = self.download_and_preprocess()
        print("Finished preprocessing dataset")

        super(AnnDataset, self).__init__(*GeneExpressionDataset.get_attributes_from_matrix(data),
                                         gene_names=gene_names)

    def preprocess(self):
        print("Preprocessing data")
        if not os.path.isfile(self.save_path + self.filename):  # nothing extracted yet
            with zipfile.ZipFile(self.save_path + self.download_name) as zf:
                print("Extracting zip file")
                zf.extractall(path=self.save_path)

        ad = anndata.read_h5ad(self.save_path + self.filename)
        gene_names = np.array(ad.obs.index.values, dtype=str)
        data = ad.X.T  # change gene * cell to cell * gene
        select = data.sum(axis=1) > 0  # Take out cells that doesn't express any gene
        data = data[select, :]

        return data, gene_names
