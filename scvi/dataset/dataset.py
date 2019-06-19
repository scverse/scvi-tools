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

    def update_genes(self, subset_genes: np.ndarray):
        """Performs a in-place sub-sampling of genes and gene-related attributes.

        Sub-selects genes according to ``subset_genes`` sub-index.
        Consequently, modifies in-place the data ``X`` and the registered gene attributes.

        :param subset_genes: Index used for gene sub-sampling.
            Either a ``int`` array with arbitrary shape which values are the indexes of the genes to keep.
            Or boolean array used as a mask-like index.
        """
        new_nb_genes = (
            len(subset_genes)
            if subset_genes.dtype is not np.dtype("bool")
            else subset_genes.sum()
        )
        logger.info(
            "Downsampling from {nb_genes} to {new_nb_genes} genes".format(
                nb_genes=self.nb_genes, new_nb_genes=new_nb_genes
            )
        )
        # update datasets
        self.X = self.X[:, subset_genes]
        for version_name in self.dataset_versions:
            data = getattr(self, version_name)
            setattr(self, version_name, data[:, subset_genes])

        # update gene-related attributes accordingly
        for attribute_name in self.gene_attribute_names:
            attr = getattr(self, attribute_name)
            setattr(self, attribute_name, attr[subset_genes])

        # remove non-expressing cells
        mask_cells_to_keep = np.asarray(self.X.sum(axis=1) > 0)
        removed_idx = np.logical_not(mask_cells_to_keep).nonzero()[0]
        if len(self.X) != len(removed_idx):
            logger.info(
                "Cells with zero expression in all genes considered were removed, "
                "the indices of the removed cells in the expression matrix were: "
                "{idxs}".format(idxs=removed_idx)
            )
        self.update_cells(mask_cells_to_keep)

    def update_cells(self, subset_cells):
        """Performs a in-place sub-sampling of cells and cell-related attributes.

        Sub-selects cells according to ``subset_cells`` sub-index.
        Consequently, modifies in-place the data ``X``, its versions and the registered cell attributes.

        :param subset_cells: Index used for cell sub-sampling.
            Either a ``int`` array with arbitrary shape which values are the indexes of the cells to keep.
            Or boolean array used as a mask-like index.
        """
        nb_cells_old = self.nb_cells

        # update gene-related attributes accordingly
        for attribute_name in self.gene_attribute_names:
            attr = getattr(self, attribute_name)
            setattr(self, attribute_name, attr[subset_cells])

        # subsample datasets which calls library_size_batch which requires already updated attributes
        self.X = self.X[subset_cells]
        for version_name in self.dataset_versions:
            data = getattr(self, version_name)
            setattr(self, version_name, data[:, subset_cells])

        logging.info(
            "Downsampled from {nb_cells} to {new_nb_cells} cells".format(
                nb_cells=nb_cells_old, new_nb_cells=self.nb_cells
            )
        )

    def subsample_genes(
        self,
        new_n_genes: int = None,
        subset_genes: Union[List[int], List[bool], np.ndarray] = None,
    ):
        """Wrapper around ``update_genes`` allowing for manual and automatic (based on count variance) subsampling."""
        # Do nothing if asked to subsample more genes than total number of genes
        if subset_genes is None and (
            new_n_genes is False or new_n_genes >= self.nb_genes
        ):
            logger.info(
                "Not subsampling since new_n_genes is None or superior to nb_genes."
            )
            return None

        if subset_genes is None:
            std_scaler = StandardScaler(with_mean=False)
            std_scaler.fit(self.X.astype(np.float64))
            subset_genes = np.argsort(std_scaler.var_)[::-1][:new_n_genes]
        self.update_genes(subset_genes)

    def subsample_cells(self, size: float = 1.0):
        """Wrapper around ``update_cells`` allowing for automatic (based on sum of counts) subsampling."""
        new_n_cells = int(size * self.nb_genes) if type(size) is not int else size
        indices = np.argsort(np.asarray(self.X.sum(axis=1)).ravel())[::-1][:new_n_cells]
        self.update_cells(indices)

    def _get_genes_filter_mask_by_attribute(
        self,
        attribute_values_to_keep: Union[List, np.ndarray],
        attribute_name: str = "gene_names",
        return_data: bool = True,
    ):
        """Returns a mask with shape (nb_genes,) equal to ``True`` if the filtering condition is.

        Specifically, the mask is true if ``self.attribute_name`` is in ``attribute_values_to_keep``.

        :param attribute_values_to_keep: Values to accept for the filtering attribute.
        :param attribute_name: Name of the attribute to filter genes on.
        :param return_data: If True, returns the filtered data along with the mask.
        """
        if attribute_name not in self.gene_attribute_names:
            raise ValueError(
                "{name} is not a registered gene attribute".format(name=attribute_name)
            )

        attribute_values = getattr(self, attribute_name)
        subset_genes = np.isin(attribute_values, attribute_values_to_keep)

        if return_data:
            return self.X[:, subset_genes], subset_genes
        else:
            return subset_genes

    def filter_genes(
        self, values_to_keep: Union[List, np.ndarray], on: str = "gene_names"
    ):
        """Performs in-place gene filtering based on any gene attribute."""
        subset_genes = self._get_genes_filter_mask_by_attribute(
            attribute_values_to_keep=values_to_keep,
            attribute_name=on,
            return_data=False,
        )
        self.update_genes(subset_genes)

    def _get_cells_filter_mask_by_attribute(
        self,
        attribute_values_to_keep: Union[List, np.ndarray],
        attribute_name: str = "labels",
        return_data: bool = True,
    ):
        """Returns a mask with shape (nb_cells,) equal to ``True`` if the filtering condition is.

        Specifically, the mask is true if ``self.attribute_name`` is in ``attribute_values_to_keep``.

        :param attribute_values_to_keep: Values to accept for the filtering attribute.
        :param attribute_name: Name of the attribute to filter cells on.
        :param return_data: If True, returns the filtered data along with the mask.
        """
        if attribute_name not in self.cell_attribute_names:
            raise ValueError(
                "{name} is not a registered cell attribute".format(name=attribute_name)
            )

        attribute_values = getattr(self, attribute_name)
        subset_cells = np.isin(attribute_values, attribute_values_to_keep)

        if return_data:
            return self.X[subset_cells], subset_cells
        else:
            return subset_cells

    def cell_types_to_label(self, cell_types: Union[List[str], np.ndarray]):
        """Forms the list of labels corresponding to the specified ``cell_types``."""
        labels = [
            np.where(getattr(self, "cell_types") == cell_type)[0][0] for cell_type in cell_types
        ]
        return np.asarray(labels, dtype=np.int64)

    def _gene_idx(self, genes):
        if type(genes[0]) is not int:
            genes_idx = [np.where(gene == getattr(self, "gene_names"))[0][0] for gene in genes]
        else:
            genes_idx = genes
        return np.asarray(genes_idx, dtype=np.int64)

    def filter_cell_types(self, cell_types: Union[List[str], List[int], np.ndarray]):
        """Performs in-place filtering of cells by keeping cell types in ``cell_types``.

        :param cell_types: numpy array of type np.int (indices) or np.str (cell-types names)
        :return:
        """
        cell_types = np.asarray(cell_types)
        if cell_types.dtype is str:
            labels = self.cell_types_to_label(cell_types)

        elif cell_types.dtype is int:
            labels = cell_types

        else:
            raise ValueError(
                "Wrong dtype for cell_types. Should be either str or int (labels)."
            )

        subset_cells = self._get_cells_filter_mask_by_attribute(
            attribute_name="labels", attribute_values_to_keep=labels, return_data=False
        )

        self.update_cells(subset_cells)

    def merge_cell_types(
        self,
        cell_types: Union[Tuple[int, ...], Tuple[str, ...], List[int], List[str], np.ndarray],
        new_cell_type_name: str = None,
    ):
        """Merges some cell types into a new one, and changes the labels accordingly.

        :param cell_types: Cell types to merge.
        :param new_cell_type_name: Name for the new aggregate cell type.
        """
        labels_subset = self.cell_types_to_label(cell_types)
        # labels should be set not muted
        new_labels = self.labels
        new_labels[
            np.isin(new_labels, labels_subset)
        ] = self.n_labels  # Put at the end the new merged cell-type
        self.labels = new_labels
        if new_cell_type_name:
            getattr(self, "cell_types")[-1] = new_cell_type_name

    def map_cell_types(
        self,
        cell_types_dict: Dict[Union[int, str, Tuple[int, ...], Tuple[str, ...]], str],
    ):
        """Performs in-place filtering of cells using a cell type mapping.

        Cell types in the keys of ``cell_types_dict`` are merged and given the name of the associated value

        :param cell_types_dict: dictionary with tuples of cell types to merge as keys
            and new cell type names as values.
        """
        keys = [
            (key,) if type(key) is not tuple else key for key in cell_types_dict.keys()
        ]
        cell_types = [cell_type for cell_types in keys for cell_type in cell_types]
        self.filter_cell_types(cell_types)
        for cell_types, new_cell_type_name in cell_types_dict.items():
            self.merge_cell_types(cell_types, new_cell_type_name)

    def corrupt(self, rate=0.1, corruption="uniform"):
        """Forms a corrupted_X attribute containing a corrupted version of X.

        Sub-samples ``rate * self.X.shape[0] * self.X.shape[1]`` entries
        and perturbs them according to the ``corruption`` method.
        Namely:
            - "uniform" multiplies the count by a Bernouilli(0.9)
            - "binomial" replaces the count with a Binomial(count, 0.2)
        A corrupted version of ``self.X`` is stored in ``self.corrupted_X``.

        :param rate: Rate of corrupted entries.
        :param corruption: Corruption method.
        """
        if (
            corruption == "uniform"
        ):  # multiply the entry n with a Ber(0.9) random variable.
            i, j = self.X.nonzero()
            ix = np.random.choice(len(i), int(np.floor(rate * len(i))), replace=False)
            i, j = i[ix], j[ix]
            self.corrupted_X[i, j] = np.multiply(
                self.X[i, j],
                np.random.binomial(n=np.ones(len(ix), dtype=np.int32), p=0.9),
            )
        elif (
            corruption == "binomial"
        ):  # multiply the entry n with a Bin(n, 0.2) random variable.
            i, j = (k.ravel() for k in np.indices(self.X.shape))
            ix = np.random.choice(len(i), int(np.floor(rate * len(i))), replace=False)
            i, j = i[ix], j[ix]
            self.corrupted_X[i, j] = np.random.binomial(
                n=(self.X[i, j]).astype(np.int32), p=0.2
            )
        else:
            raise NotImplementedError("Unknown corruption method.")

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

    def raw_counts_properties(
        self, idx1: Union[List[int], np.ndarray], idx2: Union[List[int], np.ndarray]
    ) -> Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
        """Computes and returns some statistics on the raw counts of two sub-populations.

        :param idx1: subset of indices desccribing the first population.
        :param idx2: subset of indices desccribing the second population.

        :return: Tuple of np.ndarray containing, by pair (one for each sub-population),
            mean expression, mean of non-zero expression, mean of normalized expression.
        """
        mean1 = (self.X[idx1, :]).mean(axis=0)
        mean2 = (self.X[idx2, :]).mean(axis=0)
        nonz1 = (self.X[idx1, :] != 0).mean(axis=0)
        nonz2 = (self.X[idx2, :] != 0).mean(axis=0)
        if self.norm_X is None:
            scaling_factor = self.X.mean(axis=1)
            self.norm_X = self.X / scaling_factor.reshape(len(scaling_factor), 1)
        norm_mean1 = self.norm_X[idx1, :].mean(axis=0)
        norm_mean2 = self.norm_X[idx2, :].mean(axis=0)
        return (
            np.asarray(mean1).ravel(),
            np.asarray(mean2).ravel(),
            np.asarray(nonz1).ravel(),
            np.asarray(nonz2).ravel(),
            np.asarray(norm_mean1).ravel(),
            np.asarray(norm_mean2).ravel(),
        )

    @classmethod
    def from_batch_array(cls, X: np.ndarray) -> "GeneExpressionDataset":
        """Forms a GeneExpressionDataset from an array with shape (n_batches, nb_cells, nb_genes).

        :param X: np.ndarray with shape (n_batches, nb_cells, nb_genes).
        """
        if len(np.squeeze(X.shape)) != 3:
            raise ValueError(
                "Shape of np.squeeze(X) != 3. Use standard constructor "
                "if your dataset has shape (nb_cells, nb_genes)"
            )
        batch_indices = np.arange(X.shape[0])[:, None] * np.ones(X.shape[1])[None, :]
        batch_indices = batch_indices.reshape(-1)
        X = X.reshape(-1, X.shape[2])
        return cls(X=X, batch_indices=batch_indices)

    @classmethod
    def from_per_batch_list(
        cls,
        Xs: List[Union[sp_sparse.csr_matrix, np.ndarray]],
        list_labels: List[Union[List[int], np.ndarray]] = None,
    ) -> "GeneExpressionDataset":
        """Forms a GeneExpressionDataset from a list of batch.

        :param Xs: RNA counts in the form of a list of np.ndarray with shape (..., nb_genes)
        :param list_labels: list of cell-wise labels for each batch.
        """
        nb_genes = Xs[0].shape[1]
        if not all(X.shape[1] == nb_genes for X in Xs):
            raise ValueError("All batches must have same nb_genes")

        X = np.concatenate(Xs) if type(Xs[0]) is np.ndarray else sp_sparse.vstack(Xs)
        batch_indices = np.concatenate(
            [i * np.ones(len(batch), dtype=np.int64) for i, batch in enumerate(Xs)]
        )
        labels = np.concatenate(list_labels).astype(np.int64) if list_labels else None

        return cls(X=X, batch_indices=batch_indices, labels=labels)

    @classmethod
    def from_per_label_list(
        cls,
        Xs: List[Union[sp_sparse.csr_matrix, np.ndarray]],
        list_batches: List[Union[List[int], np.ndarray]] = None,
    ) -> "GeneExpressionDataset":
        """Forms a GeneExpressionDataset from a list of batch.

        :param Xs: RNA counts in the form of a list of np.ndarray with shape (..., nb_genes)
        :param list_batches: list of cell-wise labels for each batch.
        """
        nb_genes = Xs[0].shape[1]
        if not all(X.shape[1] == nb_genes for X in Xs):
            raise ValueError("All batches must have same nb_genes")

        X = np.concatenate(Xs) if type(Xs[0]) is np.ndarray else sp_sparse.vstack(Xs)
        labels = np.concatenate(
            [i * np.ones(len(cluster), dtype=np.int64) for i, cluster in enumerate(Xs)]
        )
        batch_indices = (
            np.concatenate(list_batches).astype(np.int64) if list_batches else None
        )

        return cls(X=X, batch_indices=batch_indices, labels=labels)

    @classmethod
    def from_datasets(
        cls,
        *gene_datasets: "GeneExpressionDataset",
        on: str = "gene_names",
        sharing_intstructions_dict: Dict[str, bool] = None,
    ) -> "GeneExpressionDataset":
        """Merges multiple datasets keeping the intersection of their genes.

        By default, merging is performed using their ``gene_names`` attribute.
        Note that datasets should all have gene_dataset.n_labels=0.
        Batch indices are generated in the same order as datasets are given.

        :param gene_datasets: a sequence of ``GeneExpressionDataset`` objects.
        :param on: attribute to select gene interesection
        :param sharing_intstructions_dict:

        :return: GeneExpressionDataset instance of the concatenated datasets
        """
        # sanity check
        if not all([hasattr(gene_dataset, on) for gene_dataset in gene_datasets]):
            raise ValueError(
                "All datasets should have the merge key 'on' as an attribute"
            )

        # set default sharing behaviour
        if sharing_intstructions_dict is None:
            sharing_intstructions_dict["batch_indices"] = "offset"

        # get insterection based on gene attribute `on` and get attribute intersection
        genes_to_keep = set.intersection(
            *[set(getattr(gene_dataset, on)) for gene_dataset in gene_datasets]
        )
        gene_attributes_to_keep = set.intersection(
            *[set(gene_dataset.gene_attribute_names) for gene_dataset in gene_datasets]
        )
        cell_attributes_to_keep = set.intersection(
            *[set(gene_dataset.cell_attribute_names) for gene_dataset in gene_datasets]
        )

        # keep gene order and attributes of the first dataset
        gene_to_keep = [
            gene for gene in getattr(gene_datasets[0], on) if gene in genes_to_keep
        ]
        logger.info("Keeping {nb_genes} genes".format(nb_genes=len(genes_to_keep)))

        # filter genes
        Xs = [dataset.filter_genes(gene_to_keep, on=on).X for dataset in gene_datasets]

        # concatenate data
        if all([type(X) is np.ndarray for X in Xs]):
            X = np.concatenate(Xs)
        # if sparse, cast all to sparse and stack
        else:
            X = sp_sparse.vstack(
                [
                    X
                    if isinstance(X, sp_sparse.csr_matrix)
                    else sp_sparse.csr_matrix(X)
                    for X in Xs
                ]
            )

        # instantiate new dataset and fill attributes
        dataset = cls(X)

        # keep gene attributes of first dataset, and keep all mappings (e.g gene types)
        for attribute_name in gene_attributes_to_keep:
            dataset.initialize_gene_attribute(
                getattr(gene_datasets[0], attribute_name), attribute_name
            )
            for gene_dataset in gene_datasets:
                mapping_names = gene_dataset.attribute_mappings[attribute_name]
                for mapping_name in mapping_names:
                    dataset.initialize_mapped_attribute(
                        attribute_name,
                        mapping_name,
                        getattr(gene_dataset, mapping_name),
                    )

        # handle cell attributes
        for attribute_name in cell_attributes_to_keep:
            instruction = sharing_intstructions_dict.get(attribute_name, "concatenate")

            mapping_names_to_keep = list(
                set.intersection(
                    *[
                        set(gene_dataset.attribute_mappings[attribute_name])
                        for gene_dataset in gene_datasets
                    ]
                )
            )

            attribute_values = []

            if instruction == "concatenate":
                mappings = defaultdict(set)
                # create new mapping
                for mapping_name in mapping_names_to_keep:
                    mappings[mapping_name] = mappings[mapping_name].union(
                        getattr(dataset, mapping_name)
                    )
                mappings = {k: list(v) for k, v in mappings.items()}

                for gene_dataset in gene_datasets:
                    local_attribute_values = np.squeeze(
                        getattr(gene_dataset, attribute_name)
                    )
                    # remap attribute according to new mapping
                    if mapping_names_to_keep:
                        ref_mapping_name = mapping_names_to_keep[0]
                        old_mapping = list(getattr(gene_dataset, ref_mapping_name))
                        new_indices = [
                            mappings[ref_mapping_name].index(v) for v in old_mapping
                        ]
                        local_attribute_values, _ = remap_categories(
                            local_attribute_values, mapping_to=new_indices
                        )
                    attribute_values.append(local_attribute_values)

            elif instruction == "offset":
                mappings = defaultdict(list)
                offset = 0
                for i, gene_dataset in enumerate(gene_datasets):
                    local_attribute_values = np.squeeze(
                        getattr(gene_dataset, attribute_name)
                    )
                    attribute_values.append(offset + local_attribute_values)
                    offset += len(np.unique(local_attribute_values))
                    for mapping_name in mapping_names_to_keep:
                        mappings[mapping_name].extend(
                            getattr(gene_dataset, mapping_name)
                        )

            else:
                raise ValueError(
                    "Unknown sharing instrunction {instruction} for attribute {name}"
                    "".format(instruction=instruction, name=attribute_name)
                )

            dataset.initialize_cell_attribute(
                np.concatenate(attribute_values), attribute_name
            )

            for mapping_name, mapping_values in mappings.items():
                dataset.initialize_mapped_attribute(
                    attribute_name, mapping_name, mapping_values
                )

        return dataset


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
