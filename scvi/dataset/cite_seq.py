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
            "ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE100nnn/GSE100866/suppl/"
            "GSE100866_%s-ADT_clr-transformed.csv.gz" % s
        )
        self.urls = [url_rna, url_adt, url_adt_clr]
        # Their centered log ratio transformation for ADT counts is different from the standard clr transformation :
        # they explain they add pseudocounts (for 0 values), but do not explicit the actual transformation,
        # which doesn't seem to be simply be adding count 1 to all entries, or only 0 entries

        self.download_name_rna = '%s_rna.csv.gz' % name
        self.download_name_adt = '%s_adt.csv.gz' % name
        self.download_name_adt_centered = '%s_adt_centered.csv.gz' % name
        self.download_names = (
            self.download_name_rna,
            self.download_name_adt,
            self.download_name_adt_centered
        )

        expression_data = self.download_and_preprocess()

        super(CiteSeqDataset, self).__init__(
            *self.get_attributes_from_matrix(expression_data)
        )

    def preprocess(self):
        print("Preprocessing data")
        self.expression = expression = pd.read_csv(self.save_path + self.download_name_rna, index_col=0,
                                                   compression='gzip').T
        self.adt = adt = pd.read_csv(self.save_path + self.download_name_adt, index_col=0, compression='gzip')
        self.adt_expression = adt.T.values
        self.protein_markers = np.array(adt.index).astype(np.str)

        self.adt_centered = adt_centered = pd.read_csv(self.save_path + self.download_name_adt_centered, index_col=0,
                                                       compression='gzip')
        self.adt_expression_clr = adt_centered.T.values
        assert (self.protein_markers == np.array(adt_centered.index).astype(np.str)).all()

        gene_symbols = np.array(expression.columns, dtype=str)

        human_filter = np.array([name.startswith('HUMAN') for name in gene_symbols], dtype=np.bool)
        print("Selecting only HUMAN genes (%d / %d)" % (human_filter.sum(), len(human_filter)))
        expression_data = expression.values[:, human_filter]
        gene_symbols = gene_symbols[human_filter]

        self.gene_symbols = np.char.upper(
            np.array([name.split('_')[-1] if '_' in name else name for name in gene_symbols], dtype=np.str)
        )

        print("Finish preprocessing data")
        return expression_data


class CbmcDataset(CiteSeqDataset):
    r""" Loads cbmc dataset.

    This dataset that includes 8,617 cord blood mononuclear cells profiled using 10x along with for each cell 13
    well-characterized mononuclear antibodies. We kept the top 600 genes by variance.

    Args:
        :save_path: Save path of raw data file. Default: ``'data/'``.

    Examples:
        >>> gene_dataset = CbmcDataset()

    """
    def __init__(self, save_path='data/citeSeq/', additional_genes=600):
        self.additional_genes = additional_genes
        super(CbmcDataset, self).__init__(name="cbmc", save_path=save_path)

        # This maintains the library statistics for just umi counts
        self._X = np.ascontiguousarray(np.concatenate([self._X, self.adt_expression], axis=1), dtype=np.float32)
        self.indexes_adt = np.arange(self._X.shape[1])[-len(self.protein_markers):]

    def preprocess(self):
        print("Preprocessing cbmc data")
        # UMIs
        self.expression = expression = pd.read_csv(self.save_path + self.download_name_rna, index_col=0,
                                                   compression='gzip').T
        gene_symbols = np.array(expression.columns, dtype=str)

        human_filter = np.array([name.startswith('HUMAN') for name in gene_symbols], dtype=np.bool)
        print("Selecting only HUMAN genes (%d / %d)" % (human_filter.sum(), len(human_filter)))
        expression_data = expression.values[:, human_filter]

        # Keep top 600 genes by variance as in scVI paper
        print('Keeping top {} genes by variance'.format(self.additional_genes))
        gene_vars = np.var(expression_data, axis=0)
        gene_inds = np.argsort(gene_vars)[-self.additional_genes:]
        expression_data = expression_data[:, gene_inds]

        # Necessary so we know which cells in adt counts to keep (could be improved)
        to_keep = np.array((expression_data.sum(axis=1) > 0)).ravel()
        expression_data = expression_data[to_keep]

        gene_symbols = gene_symbols[human_filter]
        self.gene_symbols = np.char.upper(
            np.array([name.split('_')[-1] if '_' in name else name for name in gene_symbols], dtype=np.str)
        )
        # ADTs
        self.adt = adt = pd.read_csv(self.save_path + self.download_name_adt, index_col=0, compression='gzip')
        # Remove CCR5, CCR7, and CD10 due to poor enrichments
        # as done in https://satijalab.org/seurat/multimodal_vignette.html
        self.adt = self.adt.drop(['CCR5', 'CCR7', 'CD10'], axis=0)
        self.adt = self.adt.loc[:, to_keep]
        self.adt_expression = adt.T.values
        self.protein_markers = np.array(adt.index).astype(np.str)

        print("Finish preprocessing data")
        return expression_data
