import pandas as pd
import numpy as np
import urllib.request
import os
from scipy import io
from pathlib import Path
from zipfile import ZipFile


from .dataset import GeneExpressionDataset


class PbmcDataset(GeneExpressionDataset):
    datazip_url = "https://github.com/romain-lopez/scVI-reproducibility/raw/master/additional/data.zip"

    def __init__(self, save_path='data/PBMC/'):
        self.save_path = save_path

        self.download_name = "gene_info.csv"
        self.gene_names_filename = "michael_gene_names.csv"  # from romain, manually put into PBMC folder
        self.count_filename = "count.mtx"  # from romain, manually put into PBMC folder

        expression_data, gene_names = self.download_and_preprocess()

        super(PbmcDataset, self).__init__(
            *GeneExpressionDataset.get_attributes_from_matrix(
                expression_data), gene_names=gene_names)

    def download(self):
        r = urllib.request.urlopen(self.datazip_url)
        print("Downloading data.zip")

        def readIter(f, blocksize=1000):
            while True:
                data = f.read(blocksize)
                if not data:
                    break
                yield data

        directory = str(Path(self.save_path).parent) + '/'

        if not os.path.exists(directory):
            os.makedirs(directory)
        with open(directory + 'data.zip', 'wb') as f:
            for data in readIter(r):
                f.write(data)

        with ZipFile(directory + 'data.zip', 'r') as zip:
            zip.extractall(path=directory)

    def preprocess(self):
        print("Preprocessing pbmc data")

        gene_names = np.array(pd.read_csv(self.save_path + self.gene_names_filename, index_col=0)["x"], dtype=str)
        expression = pd.DataFrame(io.mmread(self.save_path + self.count_filename).T.A, columns=gene_names)

        micro_array_result = pd.read_csv(self.save_path + self.download_name)
        expression = expression[micro_array_result["ENSG"]]

        de_expression = np.copy(expression.values)
        de_gene_names = np.array(expression.columns, dtype=str)

        return de_expression, de_gene_names
