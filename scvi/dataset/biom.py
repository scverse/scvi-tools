import logging
import os
from typing import Iterable, Union

import numpy as np
import pandas as pd

from scvi.dataset.dataset import DownloadableDataset
try:
    from biom import load_table
except:
    print('`biom-format` not installed. Run `pip install biom-format` '
          'to make warning go away.')

logger = logging.getLogger(__name__)



class BiomDataset(DownloadableDataset):
    """Loads a `.biom` file.

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
        >>> remote_url = "https://github.com/biocore/American-Gut/"
        ... "blob/master/data/HMP/HMPv35_100nt.biom?raw=true"
        >>> remote_biom_dataset = BiomDataset("HMP.biom", save_path='data/',
        ... url=remote_url)
        >>> # Loading a local dataset
        >>> local_csv_dataset = CsvDataset("HMP.biom", save_path="data/")
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
            gene_by_cell  # Whether the original dataset is genes by cells
        )
        self.labels_file = labels_file
        self.batch_ids_file = batch_ids_file
        super().__init__(
            urls=url,
            filenames=filename,
            save_path=save_path,
            delayed_populating=delayed_populating,
        )
        if (new_n_genes is not None) or (subset_genes is not None):
            self.subsample_genes(new_n_genes=new_n_genes, subset_genes=subset_genes)

    def read_biom(self, path):
        """Loads a `.biom` file from a specified local path

        :param path: File name to use when loading the data.

        """
        table = load_table(path)
        # D : genes, N : samples (or cells)
        sp_matrix = table.matrix_data  # csr_matrix: D x N
        cells = table.ids(axis='sample')
        genes = table.ids(axis='observation')
        return sp_matrix, cells, genes

    def populate(self):
        logger.info("Preprocessing dataset")
        path = os.path.join(self.save_path, self.filenames[0])
        sp_matrix, _, gene_names = self.read_biom(path)

        labels, cell_types, batch_indices = None, None, None
        if self.labels_file is not None:
            labels = pd.read_csv(
                os.path.join(self.save_path, self.labels_file),
                header=0, index_col=0
            )
            labels = labels.values
            cell_types = np.unique(labels)

        if self.batch_ids_file is not None:
            batch_indices = pd.read_csv(
                os.path.join(self.save_path, self.batch_ids_file),
                header=0, index_col=0
            )
            batch_indices = batch_indices.values

        logger.info("Finished preprocessing dataset")
        print(len(gene_names), sp_matrix.shape)
        self.populate_from_data(
            X=sp_matrix.T,
            batch_indices=batch_indices,
            labels=labels,
            gene_names=gene_names,
            cell_types=cell_types,
        )
        self.filter_cells_by_count()
