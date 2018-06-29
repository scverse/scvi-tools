import pandas as pd
from zipfile import ZipFile
from .dataset import GeneExpressionDataset
from pathlib import Path


class HematoDataset(GeneExpressionDataset):
    def __init__(self, save_path='data/HEMATO/'):
        self.save_path = save_path

        # file and datazip file
        self.urls = ["https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM2388072&format=file&"
                     "file=GSM2388072%5Fbasal%5Fbone%5Fmarrow%2Eraw%5Fumifm%5Fcounts%2Ecsv%2Egz",
                     "https://github.com/romain-lopez/scVI-reproducibility/raw/master/additional/data.zip"]

        # Originally: GSM2388072_basal_bone_marrow.raw_umifm_counts.csv.gz
        self.download_names = ["bBM.raw_umifm_counts.csv.gz", "data.zip"]
        self.gene_names_filename = "bBM.filtered_gene_list.paper.txt"
        self.spring_and_pba_filename = "bBM.spring_and_pba.csv"

        expression_data, gene_names = self.download_and_preprocess()

        super(HematoDataset, self).__init__(
            *GeneExpressionDataset.get_attributes_from_matrix(
                expression_data), gene_names=gene_names)

    def preprocess(self):
        print("Preprocessing Hemato data")

        with ZipFile(self.save_path + 'data.zip', 'r') as zip:
            zip.extractall(path=Path(self.save_path).parent)

        raw_counts = pd.read_csv(self.save_path + self.download_names[0], compression='gzip')
        raw_counts.drop(raw_counts.index[raw_counts["library_id"] == "basal_bm1"], inplace=True)

        spring_and_pba = pd.read_csv(self.save_path + self.spring_and_pba_filename)
        with open(self.save_path + self.gene_names_filename) as f:
            gene_filter_list = f.read()

        gene_names = gene_filter_list.splitlines()

        data = raw_counts.merge(spring_and_pba, how="inner")
        expression_data = data[gene_names]

        expression_data = expression_data.values

        print("Finished preprocessing Hemato data")
        return expression_data, gene_names
