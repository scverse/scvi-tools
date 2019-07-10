import logging
import os
from typing import Iterable, Union

import numpy as np
import pandas as pd

from scvi.dataset.dataset import DownloadableDataset

logger = logging.getLogger(__name__)


class CsvDataset(DownloadableDataset):
    """Loads a `.csv` file.

    :param filename: File name to use when saving/loading the data.
    :param save_path: Location to use when saving/loading the data.
    :param url: URL pointing to the data which will be downloaded
        if it's not already in ``save_path``.
    :param new_n_genes: Number of subsampled genes.
    :param subset_genes: List of genes for subsampling.
    :param compression: For on-the-fly decompression of on-disk data. If ‘infer’ and filepath_or_buffer
        is path-like, then detect compression from the following extensions: ‘.gz’, ‘.bz2’, ‘.zip’, or ‘.xz’
        (otherwise no decompression). If using ‘zip’, the ZIP file must contain only one data file to be read in.
    :param batch_ids_file: Name of the `.csv` file with batch indices.
        File contains two columns. The first holds cell names and second
        holds batch indices - type int. The first row of the file is header.

    Examples:
        >>> # Loading a remote dataset
        >>> remote_url = "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE100866&format=file&file="
        ... "GSE100866%5FCBMC%5F8K%5F13AB%5F10X%2DRNA%5Fumi%2Ecsv%2Egz")
        >>> remote_csv_dataset = CsvDataset("GSE100866_CBMC_8K_13AB_10X-RNA_umi.csv.gz", save_path='data/',
        ... compression="gzip", url=remote_url)
        >>> # Loading a local dataset
        >>> local_csv_dataset = CsvDataset("GSE100866_CBMC_8K_13AB_10X-RNA_umi.csv.gz",
        ... save_path="data/", compression='gzip')
    """

    def __init__(
        self,
        filename: str,
        save_path: str = "data/",
        url: str = None,
        new_n_genes: int = None,
        subset_genes: Iterable[Union[int, str]] = None,
        compression: str = None,
        sep: str = ",",
        gene_by_cell: bool = True,
        labels_file: str = None,
        batch_ids_file: str = None,
        delayed_populating: bool = False,
    ):
        self.compression = compression
        self.sep = sep
        self.gene_by_cell = (
            gene_by_cell
        )  # Whether the original dataset is genes by cells
        self.labels_file = labels_file
        self.batch_ids_file = batch_ids_file
        super().__init__(
            urls=url,
            filenames=filename,
            save_path=save_path,
            delayed_populating=delayed_populating,
        )
        self.subsample_genes(new_n_genes, subset_genes)

    def populate(self):
        logger.info("Preprocessing dataset")

        if self.gene_by_cell:
            data = pd.read_csv(
                os.path.join(self.save_path, self.filenames[0]),
                sep=self.sep,
                index_col=0,
                compression=self.compression,
            ).T
        else:
            data = pd.read_csv(
                os.path.join(self.save_path, self.filenames[0]),
                sep=self.sep,
                index_col=0,
                compression=self.compression,
            )

        gene_names = np.asarray(data.columns, dtype=str)
        labels, cell_types, batch_indices = None, None, None
        if self.labels_file is not None:
            labels = pd.read_csv(
                os.path.join(self.save_path, self.labels_file), header=0, index_col=0
            )
            labels = labels.values
            cell_types = np.unique(labels)

        if self.batch_ids_file is not None:
            batch_indices = pd.read_csv(
                os.path.join(self.save_path, self.batch_ids_file), header=0, index_col=0
            )
            batch_indices = batch_indices.values

        data = data.values
        logger.info("Finished preprocessing dataset")

        self.populate_from_data(
            X=data,
            batch_indices=batch_indices,
            labels=labels,
            gene_names=gene_names,
            cell_types=cell_types,
        )
        self.filter_cells_by_count()


class BreastCancerDataset(CsvDataset):
    def __init__(self, save_path: str = "data/", delayed_populating: bool = False):
        super().__init__(
            "Layer2_BC_count_matrix-1.tsv",
            save_path=save_path,
            url="http://www.spatialtranscriptomicsresearch.org/wp-content/"
            "uploads/2016/07/Layer2_BC_count_matrix-1.tsv",
            sep="\t",
            gene_by_cell=False,
            delayed_populating=delayed_populating,
        )


class MouseOBDataset(CsvDataset):
    def __init__(self, save_path: str = "data/", delayed_populating: bool = False):
        super().__init__(
            "Rep11_MOB_count_matrix-1.tsv",
            save_path=save_path,
            url="http://www.spatialtranscriptomicsresearch.org/wp-content/uploads/"
            "2016/07/Rep11_MOB_count_matrix-1.tsv",
            sep="\t",
            gene_by_cell=False,
            delayed_populating=delayed_populating,
        )
