import numpy as np
import pandas as pd

from .dataset import GeneExpressionDataset

available_datasets = {
    "cbmc": "CBMC_8K_13AB",
    "pbmc": "PBMC_vs_flow"
}


class CiteSeqDataset(GeneExpressionDataset):
    def __init__(self, name='cbmc', save_path='data/citeSeq/'):
        self.save_path = save_path + name + '/'

        # originally: GSE100866_CBMC_8K_13AB_10X-RNA_umi.csv.gz
        s = available_datasets[name]
        url_rna, url_adt = (
            "ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE100nnn/GSE100866/suppl/GSE100866_%s_10X-RNA_umi.csv.gz" % s,
            "ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE100nnn/GSE100866/suppl/GSE100866_%s_10X-ADT_umi.csv.gz" % s
        )

        self.download_name_rna = '%s_rna.csv.gz' % name
        self.download_name_adt = '%s_adt.csv.gz' % name
        GeneExpressionDataset._download(url_adt, self.save_path, self.download_name_adt)
        GeneExpressionDataset._download(url_rna, self.save_path, self.download_name_rna)

        expression_data, gene_names = self.preprocess()

        super(CiteSeqDataset, self).__init__(
            *GeneExpressionDataset.get_attributes_from_matrix(
                expression_data), gene_names=gene_names)

    def preprocess(self):
        print("Preprocessing cbmc data")
        expression = pd.read_csv(self.save_path + self.download_name_rna, index_col=0, compression='gzip').T
        adt = pd.read_csv(self.save_path + self.download_name_adt, index_col=0,
                          compression='gzip')
        self.adt_expression = adt.T.values
        self.protein_markers = np.array(adt.index).astype(np.str)
        gene_names = np.array(expression.columns, dtype=str)

        selected = np.std(expression.values, axis=0).argsort()[-600:][::-1]
        expression_data = expression.values[:, selected]

        return expression_data, gene_names
