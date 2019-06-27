import logging
import os

import loompy
import numpy as np

from scvi.dataset.dataset import DownloadableDataset

logger = logging.getLogger(__name__)


class LoomDataset(DownloadableDataset):
    r""" Loads a `.loom` file.

    Examples:
        >>> # Loading a remote dataset
        >>> remote_loom_dataset = LoomDataset("osmFISH_SScortex_mouse_all_cell.loom", save_path='data/',
        ... url='http://linnarssonlab.org/osmFISH/osmFISH_SScortex_mouse_all_cells.loom')
        >>> # Loading a local dataset
        >>> local_loom_dataset = LoomDataset("osmFISH_SScortex_mouse_all_cell.loom", save_path='data/')
    """

    def __init__(self, filename: str, save_path: str = 'data/', url: str = None, delayed_populating: bool = False):
        super().__init__(
            urls=url,
            filenames=filename,
            save_path=save_path,
            delayed_populating=delayed_populating,
        )

    def populate(self):
        logger.info("Preprocessing dataset")
        gene_names, labels, batch_indices, cell_types = None, None, None, None
        ds = loompy.connect(os.path.join(self.save_path, self.filenames[0]))
        select = ds[:, :].sum(axis=0) > 0  # Take out cells that doesn't express any gene

        if 'Gene' in ds.ra:
            gene_names = ds.ra['Gene']

        if 'BatchID' in ds.ca:
            batch_indices = ds.ca['BatchID']
            batch_indices = np.reshape(batch_indices, (batch_indices.shape[0], 1))[select]

        if 'ClusterID' in ds.ca:
            labels = np.array(ds.ca['ClusterID'])
            labels = np.reshape(labels, (labels.shape[0], 1))[select]

        if 'CellTypes' in ds.attrs:
            cell_types = np.array(ds.attrs['CellTypes'])

        data = ds[:, select].T  # change matrix to cells by genes
        ds.close()

        logger.info("Finished preprocessing dataset")
        self.populate_from_data(
            X=data,
            batch_indices=batch_indices,
            labels=labels,
            gene_names=gene_names,
            cell_types=cell_types,
        )


class RetinaDataset(LoomDataset):
    r""" Loads retina dataset.

    The dataset of bipolar cells contains after their original pipeline for filtering 27,499 cells and
    13,166 genes coming from two batches. We use the cluster annotation from 15 cell-types from the author.
    We also extract their normalized data with Combat and use it for benchmarking.

    Examples:
        >>> gene_dataset = RetinaDataset()
    """
    def __init__(self, save_path='data/'):
        super().__init__(filename='retina.loom',
                         save_path=save_path,
                         url='https://github.com/YosefLab/scVI-data/raw/master/retina.loom')
        self.cell_types = ["RBC", "MG", "BC5A", "BC7", "BC6", "BC5C", "BC1A", "BC3B", "BC1B", "BC2", "BC5D", "BC3A",
                           "BC5B", "BC4", "BC8_9"]
