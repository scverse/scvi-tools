import numpy as np
import pandas as pd

from .dataset import GeneExpressionDataset

available_datasets = {
    "cbmc": "CBMC_8K_13AB_10X",
    "pbmc": "PBMC_vs_flow_10X",
    "cd8": "CD8_merged"
}


class CiteSeqDataset(GeneExpressionDataset):
    def __init__(self, name='cbmc', save_path='data/citeSeq/'):
        self.save_path = save_path + name + '/'

        s = available_datasets[name]
        url_rna, url_adt, url_adt_clr = (
            "ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE100nnn/GSE100866/suppl/GSE100866_%s-RNA_umi.csv.gz" % s,
            "ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE100nnn/GSE100866/suppl/GSE100866_%s-ADT_umi.csv.gz" % s,
            "ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE100nnn/GSE100866/suppl/GSE100866_%s-ADT_clr-transformed.csv.gz" % s
        )
        # Their centered log ratio transformation for ADT counts is different from the standard clr transformation :
        # they explain they add pseudocounts (for 0 values), but do not explicit the actual transformation,
        # which doesn't seem to be simply be adding count 1 to all entries, or only 0 entries

        self.download_name_rna = '%s_rna.csv.gz' % name
        self.download_name_adt = '%s_adt.csv.gz' % name
        self.download_name_adt_centered = '%s_adt_centered.csv.gz' % name
        GeneExpressionDataset._download(url_adt, self.save_path, self.download_name_adt)
        GeneExpressionDataset._download(url_rna, self.save_path, self.download_name_rna)
        GeneExpressionDataset._download(url_adt_clr, self.save_path, self.download_name_adt_centered)

        expression_data, gene_names = self.preprocess()

        super(CiteSeqDataset, self).__init__(
            *GeneExpressionDataset.get_attributes_from_matrix(expression_data), gene_names=gene_names
        )

    def preprocess(self):
        print("Preprocessing citeSeq data")
        expression = pd.read_csv(self.save_path + self.download_name_rna, index_col=0, compression='gzip').T
        adt = pd.read_csv(self.save_path + self.download_name_adt, index_col=0,
                          compression='gzip')
        self.adt_expression = adt.T.values
        self.protein_markers = np.array(adt.index).astype(np.str)

        adt_centered = pd.read_csv(self.save_path + self.download_name_adt_centered, index_col=0,
                                   compression='gzip')
        self.adt_expression_clr = adt_centered.T.values
        assert (adt_centered.protein_markers == np.array(adt_centered.index).astype(np.str)).all()

        gene_names = np.array(expression.columns, dtype=str)

        if any([name.startswith('HUMAN') for name in gene_names]):  # and False:
            human_filter = np.array([name.startswith('HUMAN') for name in gene_names], dtype=np.bool)
            print("Selecting %d out of %d genes" % (human_filter.sum(), len(human_filter)))
            expression_data = expression.values[:, human_filter]
            gene_names = gene_names[human_filter]
        else:
            print("No human filter")
            expression_data = expression.values

        # gene_names = gene_names[selected]
        # These are not really gene_names
        self.gene_symbols = np.char.upper(
            np.array([name.split('_')[-1] if '_' in name else name for name in gene_names], dtype=np.str))

        return expression_data, gene_names
