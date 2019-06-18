# -*- coding: utf-8 -*-

"""Handling datasets.
For the moment, is initialized with a torch Tensor of size (n_cells, nb_genes)"""
import copy
import logging
import os
import urllib.request

from abc import abstractmethod
from collections import defaultdict
from typing import Dict, List, Tuple, Union

import numpy as np
import scipy.sparse as sp_sparse
import torch
from sklearn.preprocessing import StandardScaler
from torch.utils.data import Dataset

logger = logging.getLogger(__name__)


class GeneExpressionDataset(Dataset):
    """Generic class representing RNA counts and annotation information.

    This class is scVI's base dataset class. It gives access to several
    standard attributes: counts, number of cells, number of genes, etc.
    More importantly, it implements gene-based and cell-based filtering methods.
    It also allows the storage of cell and gene annotation information,
    as well as mappings from these annotation attributes to unique identifiers.
    In order to propagate the filtering behaviour correctly through the relevant
    attributes, they are kept in registries (cell, gene, mappings) which are
    iterated through upon any filtering operation.


    :param X: RNA counts matrix, sparse format supported (e.g ``scipy.sparse.csr_matrix``).
    :param local_means: ``np.ndarray`` with shape (nb_cells,). Mean counts per batch.
        If ``None``, they are computed automatically according to ``batch_indices``.
    :param local_vars: ``np.ndarray`` with shape (nb_cells,). Variance of counts per batch.
        If ``None``, they are computed automatically according to ``batch_indices``.
    :param batch_indices: ``np.ndarray`` with shape (nb_cells,). Maps each cell to the batch
        it originates from. Note that a batch most likely refers to a specific piece
        of tissue or a specific experimental protocol.
    :param labels: ``np.ndarray`` with shape (nb_cells,). Cell-wise labels. Can be mapped
        to cell types using attribute mappings.
    :param gene_names: ``List`` or ``np.ndarray`` with length/shape (nb_genes,).
        Maps each gene to its name.
    :param cell_types: Maps each integer label in ``labels`` to a cell type.
    :param x_coord: ``np.ndarray`` with shape (nb_cells,). x-axis coordinate of each cell.
        Useful for spatial data, e.g originating a FISH-like protocol.
    :param y_coord: ``np.ndarray`` with shape (nb_cells,). y-axis coordinate of each cell.
        Useful for spatial data, e.g originating a FISH-like protocol.
    """

    def __init__(
        self,
        X: Union[np.ndarray, sp_sparse.csr_matrix],
        batch_indices: Union[List[int], np.ndarray, sp_sparse.csr_matrix] = None,
        labels: Union[List[int], np.ndarray, sp_sparse.csr_matrix] = None,
        local_means: Union[List[float], np.ndarray] = None,
        local_vars: Union[List[float], np.ndarray] = None,
        gene_names: Union[List[str], np.ndarray, sp_sparse.csr_matrix] = None,
        cell_types: Union[List[int], np.ndarray] = None,
        x_coord: Union[List[float], np.ndarray, sp_sparse.csr_matrix] = None,
        y_coord: Union[List[float], np.ndarray, sp_sparse.csr_matrix] = None,
    ):
        # registers
        self.dataset_versions = set()
        self.gene_attribute_names = set()
        self.cell_attribute_names = set()
        self.attribute_mappings = defaultdict(list)

        # set the data hidden attribute
        self._X = (
            np.ascontiguousarray(X, dtype=np.float32)
            if isinstance(X, np.ndarray)
            else X
        )

        # handle attributes with defaults
        batch_indices = np.asarray(batch_indices) if batch_indices else np.zeros(len(X))
        self._batch_indices = batch_indices
        self.cell_attribute_names.add("batch_indices")
        self.labels = np.asarray(labels) if labels else np.zeros(len(X))
        self.cell_attribute_names.add("labels")

        # handle library size, computing requires batch_indices
        if local_means is not None and local_vars is not None:
            self.initialize_cell_attribute(np.asarray(local_means), "local_means")
            self.initialize_cell_attribute(np.asarray(local_vars), "local_vars")
        elif local_means is None and local_vars is None:
            self.library_size_batch()
            self.cell_attribute_names.update(["local_means", "local_vars"])
        else:
            raise ValueError(
                "When using custom library sizes, both local_means "
                "and local_vars should be provided."
            )

        # handle optional attributes
        if x_coord is not None:
            self.initialize_cell_attribute(np.asarray(x_coord), "x_coord")
        if y_coord is not None:
            self.initialize_cell_attribute(np.asarray(y_coord), "y_coord")
        if gene_names is not None:
            self.initialize_gene_attribute(
                np.asarray(gene_names, dtype=np.str), "gene_names"
            )
        if cell_types is not None:
            self.initialize_mapped_attribute(
                "labels", "cell_types", np.asarray(cell_types, dtype=np.str)
            )

    def __len__(self):
        return self.X.shape[0]

    def __getitem__(self, idx):
        """Implements @abstractcmethod in torch.utils.data.dataset.Dataset"""
        return idx

    @property
    def X(self):
        return self._X

    @X.setter
    def X(self, X: Union[np.ndarray, sp_sparse.csr_matrix]):
        """Recomputes the library size"""
        self._X = X
        self.library_size_batch()

    @property
    def nb_cells(self) -> int:
        return self.X.shape[0]

    @property
    def nb_genes(self) -> int:
        return self.X.shape[1]

    @property
    def batch_indices(self) -> np.ndarray:
        return self._batch_indices

    @batch_indices.setter
    def batch_indices(self, batch_indices: Union[List[int], np.ndarray]):
        """Remaps batch_indices to [0, N]"""
        batch_indices = np.asarray(batch_indices, dtype=np.int64)
        batch_indices, n_batches = remap_categories(batch_indices)
        self._n_batches = n_batches
        self._batch_indices = batch_indices

    @property
    def n_batches(self) -> int:
        try:
            return self._n_batches
        except AttributeError:
            if not hasattr(self, "_batch_indices"):
                logging.error("batch_indices attribute was not set")
            raise

    @property
    def labels(self) -> np.ndarray:
        return self._labels

    @labels.setter
    def labels(self, labels: Union[List[int], np.ndarray]):
        """Ensures that labels are always mapped to [0, 1, .., n_labels] and tracks cell_types accordingly."""
        labels = np.asarray(labels, dtype=np.int64)
        # remap to [0, 1, .. , n_labels]
        new_labels, new_n_labels = remap_categories(labels)
        # keep track of cell_types
        if hasattr(self, "cell_types"):
            self.remap_cell_types(labels)

        self._labels = new_labels
        self._n_labels = new_n_labels

    def remap_cell_types(self, labels):
        """Remaps cell_types using new labels."""
        new_cell_types = []
        n_unknown_cell_types = 0
        for new_label in np.unique(labels).astype(np.int64):
            if new_label < self.n_labels:
                new_cell_types.append(getattr(self, "cell_types")[new_label])
            # if new cell_type, needs to be set elsewhere, using 'unknown_c...' in the meantime
            else:
                new_cell_types.append("unknown_cell_type_" + str(n_unknown_cell_types))
                n_unknown_cell_types += 1

    @property
    def n_labels(self) -> int:
        try:
            return self._n_labels
        except AttributeError:
            if not hasattr(self, "_labels"):
                logging.error("attribute labels was not set")
            raise

    @property
    def norm_X(self) -> Union[sp_sparse.csr_matrix, np.ndarray]:
        """Forms a normalized version of X or returns it if it's already been computed."""
        if not hasattr(self, "_norm_X"):
            scaling_factor = self.X.mean(axis=1)
            self._norm_X = self.X / scaling_factor.reshape(len(scaling_factor), 1)
            self.register_dataset_version("norm_X")
        return self._norm_X

    @norm_X.setter
    def norm_X(self, norm_X: Union[sp_sparse.csr_matrix, np.ndarray]):
        self._norm_X = norm_X

    @property
    def corrupted_X(self) -> Union[sp_sparse.csr_matrix, np.ndarray]:
        """Returns the corrupted version of X or stores a copy of X in corrupted_X if it doesn't exist."""
        if not hasattr(self, "_corrupted_X"):
            logger.info("Storing a copy of X to be corrupted ")
            self._corrupted_X = copy.deepcopy(self.X)
            self.register_dataset_version("corrupted_X")
        return self._corrupted_X

    @corrupted_X.setter
    def corrupted_X(self, corrupted_X: Union[sp_sparse.csr_matrix, np.ndarray]):
        self._corrupted_X = corrupted_X

    def register_dataset_version(self, version_name):
        """Registers a version of the dataset, e.g normalized."""
        self.dataset_versions.add(version_name)

    def initialize_cell_attribute(self, attribute, attribute_name):
        """Sets and registers a cell-wise attribute, e.g annotation information."""
        if not self.nb_cells == len(attribute):
            raise ValueError("Number of cells and length of cell attribute mismatch")
        setattr(self, attribute_name, attribute)
        self.cell_attribute_names.add(attribute_name)

    def initialize_gene_attribute(self, attribute, attribute_name):
        """Sets and registers a gene-wise attribute, e.g annotation information."""
        if not self.nb_cells == len(attribute):
            raise ValueError("Number of genes and length of gene attribute mismatch")
        setattr(self, attribute_name, attribute)
        self.gene_attribute_names.add(attribute_name)

    def initialize_mapped_attribute(
        self, source_attribute_name, mapping_name, mapping_values
    ):
        """Sets and registers and attribute mapping, e.g labels to named cell_types."""
        if not len(np.unique(getattr(self, source_attribute_name))) == len(
            mapping_values
        ):
            raise ValueError("Number of labels and cell types mismatch")
        self.attribute_mappings[source_attribute_name].append(mapping_name)
        if mapping_values:
            setattr(self, mapping_name, mapping_values)

    def library_size_batch(self):
        """Computes the library size per batch."""
        self.local_means = np.zeros((self.nb_cells, self.nb_genes))
        self.local_vars = np.zeros((self.nb_cells, self.nb_genes))
        for i_batch in range(self.n_batches):
            idx_batch = (self.batch_indices == i_batch).ravel()
            self.local_means[idx_batch], self.local_vars[idx_batch] = library_size(
                self.X[idx_batch]
            )

    def collate_fn(
        self, batch: Union[List[int], np.ndarray]
    ) -> Tuple[torch.Tensor, ...]:
        """Batch creation function to be passed to Torch's DataLoader."""
        indices = np.asarray(batch)
        X = self.X[indices]
        return self.make_tensor_batch_from_indices(X, indices)

    def collate_fn_corrupted(
        self, batch: Union[List[int], np.ndarray]
    ) -> Tuple[torch.Tensor, ...]:
        """Batch creation function to be passed to Torch's DataLoader."""
        indices = np.asarray(batch)
        X_batch = self.corrupted_X[indices]
        return self.make_tensor_batch_from_indices(X_batch, indices)

    def make_tensor_batch_from_indices(
        self, X_batch: Union[sp_sparse.csr_matrix, np.ndarray], indices: np.ndarray
    ) -> Tuple[torch.Tensor, ...]:
        """Given indices and batch X_batch, returns a full batch of ``Torch.Tensor``"""
        if isinstance(X_batch, np.ndarray):
            X_batch = torch.from_numpy(X_batch)
        else:
            X_batch = torch.from_numpy(X_batch.toarray().astype(np.float32))
        if not hasattr(self, "x_coord") or not hasattr(self, "y_coord"):
            return (
                X_batch,
                torch.from_numpy(self.local_means[indices].astype(np.float32)),
                torch.from_numpy(self.local_vars[indices].astype(np.float32)),
                torch.from_numpy(self.batch_indices[indices].astype(np.int64)),
                torch.from_numpy(self.labels[indices].astype(np.int64)),
            )
        else:
            return (
                X_batch,
                torch.from_numpy(self.local_means[indices].astype(np.float32)),
                torch.from_numpy(self.local_vars[indices].astype(np.float32)),
                torch.from_numpy(self.batch_indices[indices].astype(np.int64)),
                torch.from_numpy(self.labels[indices].astype(np.int64)),
                torch.from_numpy(getattr(self, "x_coord")[indices].astype(np.float32)),
                torch.from_numpy(getattr(self, "y_coord")[indices].astype(np.float32)),
            )

    def corrupt(self, rate=0.1, corruption="uniform"):
        '''On the fly corruption is slow, but might be optimized in pytorch. Numpy code left here.'''
        self.corrupted_X = copy.deepcopy(self.X)
        if corruption == "uniform":  # multiply the entry n with a Ber(0.9) random variable.
            i, j = np.nonzero(self.X)
            ix = np.random.choice(range(len(i)), int(np.floor(rate * len(i))), replace=False)
            i, j = i[ix], j[ix]
            corrupted = np.multiply(self.X[i, j], np.random.binomial(n=np.ones(len(ix), dtype=np.int32), p=0.9))
        elif corruption == "binomial":  # multiply the entry n with a Bin(n, 0.9) random variable.
            i, j = (k.ravel() for k in np.indices(self.X.shape))
            ix = np.random.choice(range(len(i)), int(np.floor(rate * len(i))), replace=False)
            i, j = i[ix], j[ix]
            corrupted = np.random.binomial(n=(self.X[i, j]).astype(np.int32), p=0.2)
        self.corrupted_X[i, j] = corrupted

    def update_genes(self, subset_genes):
        new_n_genes = len(subset_genes) if subset_genes.dtype is not np.dtype('bool') else subset_genes.sum()
        logging.info("Downsampling from %i to %i genes" % (self.nb_genes, new_n_genes))
        if hasattr(self, 'gene_names'):
            self.gene_names = self.gene_names[subset_genes]
        if hasattr(self, 'gene_symbols'):
            self.gene_symbols = self.gene_symbols[subset_genes]
        self._X = self.X[:, subset_genes]
        if self.norm_X is not None:
            self.norm_X = self.norm_X[:, subset_genes]
        self.nb_genes = self.X.shape[1]
        to_keep = np.array(self.X.sum(axis=1) > 0).ravel()
        if self.X.shape != self.X[to_keep].shape:
            removed_idx = []
            for i in range(len(to_keep)):
                if not to_keep[i]:
                    removed_idx.append(i)
            logging.info("Cells with zero expression in all genes considered were "
                         "removed, the indices of the removed cells "
                         "in the expression matrix were:")
            logging.info(removed_idx)
        self.update_cells(to_keep)

    def update_cells(self, subset_cells):
        new_n_cells = len(subset_cells) if subset_cells.dtype is not np.dtype('bool') else subset_cells.sum()
        logging.info("Downsampling from %i to %i cells" % (len(self), new_n_cells))
        for attr_name in [
            '_X',
            'labels',
            'batch_indices',
            'local_means',
            'local_vars',
            'x_coord',
            'y_coord'
        ]:
            if getattr(self, attr_name) is not None:
                setattr(self, attr_name, getattr(self, attr_name)[subset_cells])
        self.library_size_batch()

    def subsample_genes(self, new_n_genes=None, subset_genes=None):
        n_cells, n_genes = self.X.shape
        if subset_genes is None and (new_n_genes is False or new_n_genes >= n_genes):
            return None  # Do nothing if subsample more genes than total number of genes
        if subset_genes is None:
            std_scaler = StandardScaler(with_mean=False)
            std_scaler.fit(self.X.astype(np.float64))
            subset_genes = np.argsort(std_scaler.var_)[::-1][:new_n_genes]
        self.update_genes(subset_genes)

    def filter_genes(self, gene_names_ref, on='gene_names'):
        """
        Same as _filter_genes but overwrites on current dataset instead of returning data,
        and updates genes names and symbols
        """
        _, subset_genes = GeneExpressionDataset._filter_genes(self, gene_names_ref, on=on)
        self.update_genes(subset_genes)

    def subsample_cells(self, size=1.):
        n_cells, n_genes = self.X.shape
        new_n_cells = int(size * n_genes) if type(size) is not int else size
        indices = np.argsort(np.array(self.X.sum(axis=1)).ravel())[::-1][:new_n_cells]
        self.update_cells(indices)

    def _cell_type_idx(self, cell_types):
        if type(cell_types[0]) is not int:
            cell_types_idx = [np.where(cell_type == self.cell_types)[0][0] for cell_type in cell_types]
        else:
            cell_types_idx = cell_types
        return np.array(cell_types_idx, dtype=np.int64)

    def _gene_idx(self, genes):
        if type(genes[0]) is not int:
            genes_idx = [np.where(gene == self.gene_names)[0][0] for gene in genes]
        else:
            genes_idx = genes
        return np.array(genes_idx, dtype=np.int64)

    def filter_cell_types(self, cell_types):
        """
        :param cell_types: numpy array of type np.int (indices) or np.str (cell-types names)
        :return:
        """
        cell_types_idx = self._cell_type_idx(cell_types)
        if hasattr(self, 'cell_types'):
            self.cell_types = self.cell_types[cell_types_idx]
            logging.info("Only keeping cell types: \n" + '\n'.join(list(self.cell_types)))
        idx_to_keep = []
        for idx in cell_types_idx:
            idx_to_keep += [np.where(self.labels == idx)[0]]
        self.update_cells(np.concatenate(idx_to_keep))
        self.labels, self.n_labels = remap_categories(self.labels, mapping_from=cell_types_idx)

    def merge_cell_types(self, cell_types, new_cell_type_name):
        """
        Merge some cell types into a new one, a change the labels accordingly.
        :param merge_cell_types: numpy array of type np.int (indices) or np.str (cell-types names)
        :return:
        """
        cell_types_idx = self._cell_type_idx(cell_types)
        for idx_from in zip(cell_types_idx):
            self.labels[self.labels == idx_from] = len(self.labels)  # Put at the end the new merged cell-type
        self.labels, self.n_labels = remap_categories(self.labels)
        if hasattr(self, 'cell_types') and type(cell_types[0]) is not int:
            new_cell_types = list(self.cell_types)
            for cell_type in cell_types:
                new_cell_types.remove(cell_type)
            new_cell_types.append(new_cell_type_name)
            self.cell_types = np.array(new_cell_types)

    def map_cell_types(self, cell_types_dict):
        """
        A map for the cell types to keep, and optionally merge together under a new name (value in the dict)
        :param cell_types_dict: a dictionary with tuples (str or int) as input and value (str or int) as output
        """
        keys = [(key,) if type(key) is not tuple else key for key in cell_types_dict.keys()]
        cell_types = [cell_type for cell_types in keys for cell_type in cell_types]
        self.filter_cell_types(cell_types)
        for cell_types, new_cell_type_name in cell_types_dict.items():
            self.merge_cell_types(cell_types, new_cell_type_name)

    def download(self):
        if hasattr(self, 'urls') and hasattr(self, 'download_names'):
            for url, download_name in zip(self.urls, self.download_names):
                GeneExpressionDataset._download(url, self.save_path, download_name)
        elif hasattr(self, 'url') and hasattr(self, 'download_name'):
            GeneExpressionDataset._download(self.url, self.save_path, self.download_name)

    @staticmethod
    def _download(url, save_path, download_name):
        if os.path.exists(os.path.join(save_path, download_name)):
            logging.info("File %s already downloaded" % (os.path.join(save_path, download_name)))
            return
        if url is None:
            logging.info("You are trying to load a local file named %s and located at %s "
                         "but this file was not found at the location %s" % (download_name, save_path,
                                                                             os.path.join(save_path, download_name)))
        r = urllib.request.urlopen(url)
        logging.info("Downloading file at %s" % os.path.join(save_path, download_name))

        def readIter(f, blocksize=1000):
            """Given a file 'f', returns an iterator that returns bytes of
            size 'blocksize' from the file, using read()."""
            while True:
                data = f.read(blocksize)
                if not data:
                    break
                yield data

        # Create the path to save the data
        if not os.path.exists(save_path):
            os.makedirs(save_path)

        with open(os.path.join(save_path, download_name), 'wb') as f:
            for data in readIter(r):
                f.write(data)

    def library_size_batch(self):
        for i_batch in range(self.n_batches):
            idx_batch = (self.batch_indices == i_batch).ravel()
            self.local_means[idx_batch], self.local_vars[idx_batch] = self.library_size(self.X[idx_batch])

    def raw_counts_properties(self, idx1, idx2):
        mean1 = (self.X[idx1, :]).mean(axis=0)
        mean2 = (self.X[idx2, :]).mean(axis=0)
        nonz1 = (self.X[idx1, :] != 0).mean(axis=0)
        nonz2 = (self.X[idx2, :] != 0).mean(axis=0)
        if self.norm_X is None:
            scaling_factor = self.X.mean(axis=1)
            self.norm_X = self.X / scaling_factor.reshape(len(scaling_factor), 1)
        norm_mean1 = self.norm_X[idx1, :].mean(axis=0)
        norm_mean2 = self.norm_X[idx2, :].mean(axis=0)
        return np.asarray(mean1).ravel(), np.asarray(mean2).ravel(), np.asarray(nonz1).ravel(), \
            np.asarray(nonz2).ravel(), np.asarray(norm_mean1).ravel(), np.asarray(norm_mean2).ravel()

    @staticmethod
    def library_size(X):
        log_counts = np.log(X.sum(axis=1))
        local_mean = (np.mean(log_counts) * np.ones((X.shape[0], 1))).astype(np.float32)
        local_var = (np.var(log_counts) * np.ones((X.shape[0], 1))).astype(np.float32)
        return local_mean, local_var

    @staticmethod
    def get_attributes_from_matrix(X, batch_indices=0, labels=None):
        ne_cells = X.sum(axis=1) > 0
        to_keep = np.where(ne_cells)[0]
        if not ne_cells.all():
            X = X[to_keep]
            removed_idx = np.where(~ne_cells)[0]
            logging.info("Cells with zero expression in all genes considered were removed, "
                         "the indices of the removed cells in the expression matrix were:")
            logging.info(removed_idx)
        local_mean, local_var = GeneExpressionDataset.library_size(X)
        batch_indices = batch_indices * np.ones((X.shape[0], 1)) if type(batch_indices) is int \
            else batch_indices[to_keep]
        labels = labels[to_keep].reshape(-1, 1) if labels is not None else np.zeros_like(batch_indices)
        return X, local_mean, local_var, batch_indices, labels

    @staticmethod
    def get_attributes_from_list(Xs, list_batches=None, list_labels=None):
        nb_genes = Xs[0].shape[1]
        assert all(X.shape[1] == nb_genes for X in Xs), "All tensors must have same size"

        new_Xs = []
        local_means = []
        local_vars = []
        batch_indices = []
        labels = []
        for i, X in enumerate(Xs):
            to_keep = np.array((X.sum(axis=1) > 0)).ravel()
            if X.shape != X[to_keep].shape:
                removed_idx = []
                for i in range(len(to_keep)):
                    if not to_keep[i]:
                        removed_idx.append(i)
                logging.info(
                    "Cells with zero expression in all genes considered were removed, the indices of the removed "
                    "cells in the ", i, "th expression matrix were:")
                logging.info(removed_idx)
            X = X[to_keep]
            new_Xs += [X]
            local_mean, local_var = GeneExpressionDataset.library_size(X)
            local_means += [local_mean]
            local_vars += [local_var]
            batch_indices += [list_batches[i][to_keep] if list_batches is not None else i * np.ones((X.shape[0], 1))]
            labels += [list_labels[i][to_keep] if list_labels is not None else np.zeros((X.shape[0], 1))]

        X = np.concatenate(new_Xs) if type(new_Xs[0]) is np.ndarray else sp_sparse.vstack(new_Xs)
        batch_indices = np.concatenate(batch_indices)
        local_means = np.concatenate(local_means)
        local_vars = np.concatenate(local_vars)
        labels = np.concatenate(labels)
        return X, local_means, local_vars, batch_indices, labels

    @staticmethod
    def concat_datasets(*gene_datasets, on='gene_names', shared_labels=True, shared_batches=False):
        """
        Combines multiple unlabelled gene_datasets based on the intersection of gene names intersection.
        Datasets should all have gene_dataset.n_labels=0.
        Batch indices are generated in the same order as datasets are given.
        :param gene_datasets: a sequence of gene_datasets object
        :return: a GeneExpressionDataset instance of the concatenated datasets
        """
        assert all([hasattr(gene_dataset, on) for gene_dataset in gene_datasets])

        gene_names_ref = set.intersection(*[set(getattr(gene_dataset, on)) for gene_dataset in gene_datasets])
        # keep gene order of the first dataset
        gene_names_ref = [gene_name for gene_name in getattr(gene_datasets[0], on) if gene_name in gene_names_ref]
        logging.info("Keeping %d genes" % len(gene_names_ref))

        Xs = [GeneExpressionDataset._filter_genes(dataset, gene_names_ref, on=on)[0] for dataset in gene_datasets]
        if gene_datasets[0].dense:
            X = np.concatenate([X if type(X) is np.ndarray else X.A for X in Xs])
        else:
            X = sp_sparse.vstack([X if type(X) is not np.ndarray else sp_sparse.csr_matrix(X) for X in Xs])

        batch_indices = np.zeros((X.shape[0], 1))
        n_batch_offset = 0
        current_index = 0
        for gene_dataset in gene_datasets:
            next_index = current_index + len(gene_dataset)
            batch_indices[current_index:next_index] = gene_dataset.batch_indices + n_batch_offset
            n_batch_offset += (gene_dataset.n_batches if not shared_batches else 0)
            current_index = next_index

        cell_types = None
        if shared_labels:
            if all([hasattr(gene_dataset, "cell_types") for gene_dataset in gene_datasets]):
                cell_types = list(
                    set([cell_type for gene_dataset in gene_datasets for cell_type in gene_dataset.cell_types])
                )
                labels = []
                for gene_dataset in gene_datasets:
                    mapping = [cell_types.index(cell_type) for cell_type in gene_dataset.cell_types]
                    labels += [remap_categories(gene_dataset.labels, mapping_to=mapping)[0]]
                labels = np.concatenate(labels)
            else:
                labels = np.concatenate([gene_dataset.labels for gene_dataset in gene_datasets])
        else:
            labels = np.zeros((X.shape[0], 1))
            n_labels_offset = 0
            current_index = 0
            for gene_dataset in gene_datasets:
                next_index = current_index + len(gene_dataset)
                labels[current_index:next_index] = gene_dataset.labels + n_labels_offset
                n_labels_offset += gene_dataset.n_labels
                current_index = next_index

        local_means = np.concatenate([gene_dataset.local_means for gene_dataset in gene_datasets])
        local_vars = np.concatenate([gene_dataset.local_vars for gene_dataset in gene_datasets])
        result = GeneExpressionDataset(X, local_means, local_vars, batch_indices, labels,
                                       gene_names=gene_names_ref, cell_types=cell_types)
        result.barcodes = [gene_dataset.barcodes if hasattr(gene_dataset, 'barcodes') else None
                           for gene_dataset in gene_datasets]
        return result

    @staticmethod
    def _filter_genes(gene_dataset, gene_names_ref, on='gene_names'):
        """
        :return: gene_dataset.X filtered by the corresponding genes ( / columns / features), idx_genes
        """
        gene_names = list(getattr(gene_dataset, on))
        subset_genes = np.array([gene_names.index(gene_name) for gene_name in gene_names_ref], dtype=np.int64)
        return gene_dataset.X[:, subset_genes], subset_genes


def remap_categories(
    original_categories: Union[List, np.ndarray],
    mapping_from: Union[List[int], np.ndarray] = None,
    mapping_to: Union[List[int], np.ndarray] = None,
) -> Tuple[np.ndarray, int]:
    """Performs and returns a remapping of a categorical array.

    If None, ``mapping_from`` is set to np.unique(categories).
    If None, ``mapping_to`` is set to [0, ..., N-1] where N is the number of categories.
    Then, ``mapping_from`` is mapped to ``mapping_to``.


    :param original_categories: Categorical array to be remapped.
    :param mapping_from: source of the mapping.
    :param mapping_to: destination of the mapping

    :return: ``tuple`` of a ``np.ndarray`` containing the new categories
        and an ``int`` equal to the new number of categories.
    """
    original_categories = np.asarray(original_categories)
    unique_categories = np.unique(original_categories)
    n_categories = len(unique_categories)
    if mapping_to is None:
        mapping_to = range(n_categories)
    if mapping_from is None:
        mapping_from = unique_categories

    # check lenghts, allow absent cell types
    if not n_categories <= len(mapping_from):
        raise ValueError(
            "Size of provided mapping_from greater than the number of categories."
        )
    if not len(mapping_to) == len(mapping_from):
        raise ValueError("Length mismatch between mapping_to and mapping_from.")

    new_categories = np.copy(original_categories)
    for idx_from, idx_to in zip(mapping_from, mapping_to):
        new_categories[original_categories == idx_from] = idx_to
    return new_categories.astype(int), n_categories


def library_size(
    X: Union[sp_sparse.csr_matrix, np.ndarray]
) -> Tuple[np.ndarray, np.ndarray]:
    log_counts = np.log(X.sum(axis=1))
    local_mean = (np.mean(log_counts) * np.ones((X.shape[0], 1))).astype(np.float32)
    local_var = (np.var(log_counts) * np.ones((X.shape[0], 1))).astype(np.float32)
    return local_mean, local_var
