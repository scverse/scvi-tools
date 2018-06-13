import pandas as pd
import os
from zipfile import ZipFile
import urllib.request
from pathlib import Path
from .dataset import GeneExpressionDataset


class HematoDataset(GeneExpressionDataset):
    url = "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM2388072&format=file&" + \
          "file=GSM2388072%5Fbasal%5Fbone%5Fmarrow%2Eraw%5Fumifm%5Fcounts%2Ecsv%2Egz"
    datazip_url = "https://github.com/romain-lopez/scVI-reproducibility/raw/master/additional/data.zip"

    def __init__(self, save_path='data/HEMATO/'):
        self.save_path = save_path

        # Originally: GSM2388072_basal_bone_marrow.raw_umifm_counts.csv.gz

        self.download_name = "bBM.raw_umifm_counts.csv.gz"
        self.gene_names_filename = "bBM.filtered_gene_list.paper.txt"
        self.spring_and_pba_filename = "bBM.spring_and_pba.csv"

        if not os.path.isfile(str(Path(save_path).parent) + '/data.zip'):
            self.download_datazip()

        expression_data, gene_names = self.download_and_preprocess()

        super(HematoDataset, self).__init__(
            *GeneExpressionDataset.get_attributes_from_matrix(
                expression_data), gene_names=gene_names)

    def download_datazip(self):
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
        print("Preprocessing Hemato data")

        raw_counts = pd.read_csv(self.save_path + self.download_name, compression='gzip')
        raw_counts.drop(raw_counts.index[raw_counts["library_id"] == "basal_bm1"], inplace=True)

        spring_and_pba = pd.read_csv(self.save_path + self.spring_and_pba_filename)
        with open(self.save_path + self.gene_names_filename) as f:
            gene_filter_list = f.read()

        gene_names = gene_filter_list.splitlines()

        data = raw_counts.merge(spring_and_pba, how="inner")
        expression_data = data[gene_names]

        expression_data = expression_data.values

        return expression_data, gene_names
