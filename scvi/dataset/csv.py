from .dataset import GeneExpressionDataset
import pandas as pd
import numpy as np
import os


class CsvDataset(GeneExpressionDataset):
    r""" Loads a `.csv` file.

    Args:
        :filename: Name of the `.csv` file.
        :save_path: Save path of the dataset. Default: ``'data/'``.
        :url: Url of the remote dataset. Default: ``None``.
        :new_n_genes: Number of subsampled genes. Default: ``600``.
        :subset_genes: List of genes for subsampling. Default: ``None``.
        :compression: For on-the-fly decompression of on-disk data. If ‘infer’ and filepath_or_buffer
            is path-like, then detect compression from the following extensions: ‘.gz’, ‘.bz2’, ‘.zip’, or ‘.xz’
            (otherwise no decompression). If using ‘zip’, the ZIP file must contain only one data file to be read in.
            Default: ``None``.
        :batch_ids_file: Name of the `.csv` file with batch indices.
            File contains two columns. The first holds gene names and second
            holds batch indices - type int. The first row of the file is header.

    Examples:
        >>> # Loading a remote dataset
        >>> remote_url = "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE100866&format=file&file=" \
        ... "GSE100866%5FCBMC%5F8K%5F13AB%5F10X%2DRNA%5Fumi%2Ecsv%2Egz")
        >>> remote_csv_dataset = CsvDataset("GSE100866_CBMC_8K_13AB_10X-RNA_umi.csv.gz", save_path='data/',
        ... compression='gzip', url=remote_url)
        >>> # Loading a local dataset
        >>> local_csv_dataset = CsvDataset("GSE100866_CBMC_8K_13AB_10X-RNA_umi.csv.gz",
        ... save_path='data/', compression='gzip')

    """

    def __init__(self, filename, save_path='data/', url=None, new_n_genes=600, subset_genes=None,
                 compression=None, sep=',', gene_by_cell=True, labels_file=None,
                 batch_ids_file=None):
        self.download_name = filename  # The given csv file is
        self.save_path = save_path
        self.url = url
        self.compression = compression
        self.sep = sep
        self.gene_by_cell = gene_by_cell  # Whether the original dataset is genes by cells
        self.labels_file = labels_file
        self.batch_ids_file = batch_ids_file

        data, gene_names, labels, cell_types, batch_ids = self.download_and_preprocess()

        super().__init__(
            *GeneExpressionDataset.get_attributes_from_matrix(
                data, labels=labels,
                batch_indices=batch_ids if batch_ids is not None else 0),
            gene_names=gene_names, cell_types=cell_types)

        self.subsample_genes(new_n_genes, subset_genes)

    def preprocess(self):
        print("Preprocessing dataset")

        if self.gene_by_cell:
            data = pd.read_csv(os.path.join(self.save_path, self.download_name),
                               sep=self.sep, index_col=0, compression=self.compression).T
        else:
            data = pd.read_csv(os.path.join(self.save_path, self.download_name),
                               sep=self.sep, index_col=0, compression=self.compression)

        gene_names = np.array(data.columns, dtype=str)
        labels, cell_types, batch_ids = None, None, None
        if self.labels_file is not None:
            labels = pd.read_csv(os.path.join(self.save_path, self.labels_file), header=0, index_col=0)
            labels = labels.values
            cell_types = np.unique(labels)

        if self.batch_ids_file is not None:
            batch_ids = pd.read_csv(
                os.path.join(
                    self.save_path, self.batch_ids_file), header=0, index_col=0)
            batch_ids = batch_ids.values

        data = data.values
        print("Finished preprocessing dataset")
        return data, gene_names, labels, cell_types, batch_ids


class BreastCancerDataset(CsvDataset):
    def __init__(self, save_path='data/'):
        super().__init__("Layer2_BC_count_matrix-1.tsv", save_path=save_path,
                         url="http://www.spatialtranscriptomicsresearch.org/wp-content/"
                             "uploads/2016/07/Layer2_BC_count_matrix-1.tsv",
                         sep='\t', gene_by_cell=False)


class MouseOBDataset(CsvDataset):
    def __init__(self, save_path='data/'):
        super().__init__("Rep11_MOB_count_matrix-1.tsv", save_path=save_path,
                         url="http://www.spatialtranscriptomicsresearch.org/wp-content/uploads/"
                             "2016/07/Rep11_MOB_count_matrix-1.tsv",
                         sep='\t', gene_by_cell=False)
