from .dataset import GeneExpressionDataset
import pandas as pd
import numpy as np


class CsvDataset(GeneExpressionDataset):
    r""" Loads a `.csv` file.

    Args:
        filename (str): Name of the `.csv` file
        save_path (str, optional): Save path of the dataset
        url (str, optional): Url of the remote dataset
        new_n_genes (int, optional): Number of subsampled genes
        subset_genes (list, optional): List of genes for subsampling
        compression(str, optional): For on-the-fly decompression of on-disk data. If ‘infer’ and filepath_or_buffer
        is path-like, then detect compression from the following extensions: ‘.gz’, ‘.bz2’, ‘.zip’, or ‘.xz’
        (otherwise no decompression). If using ‘zip’, the ZIP file must contain only one data file to be read in.

    Examples:
        >>> # Loading a remote dataset
        >>> remote_url = "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE100866&format=file&file=" \
        ... "GSE100866%5FCBMC%5F8K%5F13AB%5F10X%2DRNA%5Fumi%2Ecsv%2Egz")
        >>> remote_csv_dataset = CsvDataset("GSE100866_CBMC_8K_13AB_10X-RNA_umi.csv.gz", save_path='data/',
        ... compression='gzip', url=remove_url)
        >>> # Loading a local dataset
        >>> local_csv_dataset = CsvDataset("GSE100866_CBMC_8K_13AB_10X-RNA_umi.csv.gz",
        ... save_path='data/', compression='gzip')

    """
    def __init__(self, filename, save_path='data/', url=None, new_n_genes=600, subset_genes=None, compression=None):
        self.download_name = filename  # The given csv file is
        self.save_path = save_path
        self.url = url
        self.compression = compression

        data, gene_names = self.download_and_preprocess()

        super(CsvDataset, self).__init__(
            *GeneExpressionDataset.get_attributes_from_matrix(
                data), gene_names=gene_names)

        self.subsample_genes(new_n_genes, subset_genes)

    def preprocess(self):
        print("Preprocessing dataset")

        data = pd.read_csv(self.save_path + self.download_name, index_col=0, compression=self.compression).T
        gene_names = np.array(data.columns, dtype=str)

        data = data.values
        print("Finished preprocessing dataset")
        return data, gene_names
