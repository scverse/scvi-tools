import loompy
import numpy as np
import os
from .dataset import GeneExpressionDataset


class LoomDataset(GeneExpressionDataset):
    r""" Loads a `.loom` file.

    Args:
        :filename: Name of the `.loom` file.
        :save_path: Save path of the dataset. Default: ``'data/'``.
        :url: Url of the remote dataset. Default: ``None``.

    Examples:
        >>> # Loading a remote dataset
        >>> remote_loom_dataset = LoomDataset("osmFISH_SScortex_mouse_all_cell.loom", save_path='data/',
        ... url='http://linnarssonlab.org/osmFISH/osmFISH_SScortex_mouse_all_cells.loom')
        >>> # Loading a local dataset
        >>> local_loom_dataset = LoomDataset("osmFISH_SScortex_mouse_all_cell.loom", save_path='data/')

    """

    def __init__(self, filename, save_path='data/', url=None):
        self.download_name = filename
        self.save_path = save_path
        self.url = url

        self.has_gene, self.has_batch, self.has_cluster = False, False, False

        data, batch_indices, labels, gene_names, cell_types = self.download_and_preprocess()

        X, local_means, local_vars, batch_indices_, labels = \
            GeneExpressionDataset.get_attributes_from_matrix(data, labels=labels)
        batch_indices = batch_indices if batch_indices is not None else batch_indices_
        super().__init__(X, local_means, local_vars, batch_indices, labels,
                         gene_names=gene_names, cell_types=cell_types)

    def preprocess(self):
        print("Preprocessing dataset")
        gene_names, labels, batch_indices, cell_types = None, None, None, None
        ds = loompy.connect(os.path.join(self.save_path, self.download_name))
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

        print("Finished preprocessing dataset")
        return data, batch_indices, labels, gene_names, cell_types


class RetinaDataset(LoomDataset):
    r""" Loads retina dataset.

    The dataset of bipolar cells contains after their original pipeline for filtering 27,499 cells and
    13,166 genes coming from two batches. We use the cluster annotation from 15 cell-types from the author. We also
    extract their normalized data with Combat and use it for benchmarking.

    Args:
        :save_path: Save path of raw data file. Default: ``'data'``.

    Examples:
        >>> gene_dataset = RetinaDataset()

    """
    def __init__(self, save_path='data/'):
        super().__init__(filename='retina.loom',
                         save_path=save_path,
                         url='https://github.com/YosefLab/scVI-data/raw/master/retina.loom')
        self.cell_types = ["RBC", "MG", "BC5A", "BC7", "BC6", "BC5C", "BC1A", "BC3B", "BC1B", "BC2", "BC5D", "BC3A",
                           "BC5B", "BC4", "BC8_9"]
