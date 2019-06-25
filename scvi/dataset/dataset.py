# -*- coding: utf-8 -*-

"""Handling datasets.
For the moment, is initialized with a torch Tensor of size (n_cells, nb_genes)"""
import copy
import logging
import os
import urllib.request

from abc import abstractmethod, ABC
from collections import defaultdict
from typing import Dict, List, Tuple, Union, Optional

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

    Note that the constructor merely instantiates the GeneExpressionDataset objects.
    It should be used in combination with one of the populating method.
    Either:
        * ``populate_from_data``: to populate using a (nb_cells, nb_genes) matrix.
        * ``populate_from_per_batch_array``: to populate using a (n_batches, nb_cells, nb_genes) matrix.
        * ``populate_from_per_batch_list``: to populate using a ``n_batches``-long
            ``list`` of (nb_cells, nb_genes) matrices.
        * ``populate_from_datasets``: to populate using multiple ``GeneExperessionDataset`` objects,
            merged using the intersection of a gene-wise attribute (``gene_names`` by default).
    """

    def __init__(self):
        # registers
        self.dataset_versions = set()
        self.gene_attribute_names = set()
        self.cell_attribute_names = set()
        self.cell_categorical_attribute_names = set()
        self.attribute_mappings = defaultdict(list)

        # initialize attributes
        self._X = None
        self._batch_indices = None
        self._labels = None
        self.n_batches = None
        self.n_labels = None
        self.gene_names = None
        self.cell_types = None
        self.local_means = None
        self.local_vars = None
        self._norm_X = None
        self._corrupted_X = None

    def populate_from_data(
        self,
        X: Union[np.ndarray, sp_sparse.csr_matrix],
        batch_indices: Union[List[int], np.ndarray, sp_sparse.csr_matrix] = None,
        labels: Union[List[int], np.ndarray, sp_sparse.csr_matrix] = None,
        gene_names: Union[List[str], np.ndarray] = None,
        cell_types: Union[List[int], np.ndarray] = None,
        cell_attributes_dict: Dict[str, Union[List, np.ndarray]] = None,
        gene_attributes_dict: Dict[str, Union[List, np.ndarray]] = None,
    ):
        """Populates the data attributes of a GeneExpressionDataset object from a (nb_cells, nb_genes) matrix.

        :param X: RNA counts matrix, sparse format supported (e.g ``scipy.sparse.csr_matrix``).
        :param batch_indices: ``np.ndarray`` with shape (nb_cells,). Maps each cell to the batch
            it originates from. Note that a batch most likely refers to a specific piece
            of tissue or a specific experimental protocol.
        :param labels: ``np.ndarray`` with shape (nb_cells,). Cell-wise labels. Can be mapped
            to cell types using attribute mappings.
        :param gene_names: ``List`` or ``np.ndarray`` with length/shape (nb_genes,).
            Maps each gene to its name.
        :param cell_types: Maps each integer label in ``labels`` to a cell type.
        :param cell_attributes_dict: ``List`` or ``np.ndarray`` with shape (nb_cells,).
        :param gene_attributes_dict: ``List`` or ``np.ndarray`` with shape (nb_genes,).
        """
        # set the data hidden attribute
        self._X = (
            np.ascontiguousarray(X, dtype=np.float32)
            if isinstance(X, np.ndarray)
            else X
        )

        self.initialize_cell_attribute(
            "batch_indices",
            np.asarray(batch_indices).reshape((-1, 1))
            if batch_indices is not None
            else np.zeros((X.shape[0], 1)),
            categorical=True,
        )
        self.initialize_cell_attribute(
            "labels",
            np.asarray(labels).reshape((-1, 1))
            if labels is not None
            else np.zeros((X.shape[0], 1)),
            categorical=True,
        )

        self.compute_library_size_batch()

        if gene_names is not None:
            self.initialize_gene_attribute(
                "gene_names", np.char.lower(np.asarray(gene_names, dtype=np.str))
            )
        if cell_types is not None:
            self.initialize_mapped_attribute(
                "labels", "cell_types", np.asarray(cell_types, dtype=np.str)
            )

        # handle additional attributes
        if cell_attributes_dict:
            for attribute_name, attribute_value in cell_attributes_dict.items():
                self.initialize_cell_attribute(attribute_name, attribute_value)
        if gene_attributes_dict:
            for attribute_name, attribute_value in gene_attributes_dict.items():
                self.initialize_gene_attribute(attribute_name, attribute_value)

        self.remap_categorical_attributes()

    def populate_from_per_batch_array(
        self,
        X: np.ndarray,
        labels_per_batch: Union[np.ndarray, List[Union[List[int], np.ndarray]]] = None,
        gene_names: Union[List[str], np.ndarray] = None,
    ):
        """Populates the data attributes of a GeneExpressionDataset object
        from an array with shape (n_batches, nb_cells, nb_genes).

        :param X: np.ndarray with shape (n_batches, nb_cells, nb_genes).
        :param labels_per_batch: cell-wise labels for each batch.
        :param gene_names: gene names, stored as ``str``.
        """
        if len(np.squeeze(X.shape)) != 3:
            raise ValueError(
                "Shape of np.squeeze(X) != 3. Use populate_from_data "
                "if your dataset has shape (nb_cells, nb_genes)"
            )
        batch_indices = np.arange(X.shape[0])[:, None] * np.ones(X.shape[1])[None, :]
        batch_indices = batch_indices.reshape(-1)
        labels = (
            np.concatenate(labels_per_batch).astype(np.int64)
            if labels_per_batch is not None
            else None
        )
        X = X.reshape(-1, X.shape[2])
        self.populate_from_data(
            X=X, batch_indices=batch_indices, labels=labels, gene_names=gene_names
        )

    def populate_from_per_batch_list(
        self,
        Xs: List[Union[sp_sparse.csr_matrix, np.ndarray]],
        labels_per_batch: List[Union[List[int], np.ndarray]] = None,
        gene_names: Union[List[str], np.ndarray] = None,
    ):
        """Populates the data attributes of a GeneExpressionDataset object from a ``n_batches``-long
            ``list`` of (nb_cells, nb_genes) matrices.

        :param Xs: RNA counts in the form of a list of np.ndarray with shape (..., nb_genes)
        :param labels_per_batch: list of cell-wise labels for each batch.
        :param gene_names: gene names, stored as ``str``.
        """
        nb_genes = Xs[0].shape[1]
        if not all(X.shape[1] == nb_genes for X in Xs):
            raise ValueError("All batches must have same nb_genes")

        X = np.concatenate(Xs) if type(Xs[0]) is np.ndarray else sp_sparse.vstack(Xs)
        batch_indices = np.concatenate(
            [i * np.ones(batch.shape[0], dtype=np.int64) for i, batch in enumerate(Xs)]
        )
        labels = (
            np.concatenate(labels_per_batch).astype(np.int64)
            if labels_per_batch
            else None
        )

        self.populate_from_data(
            X=X, batch_indices=batch_indices, labels=labels, gene_names=gene_names
        )

    def populate_from_per_label_list(
        self,
        Xs: List[Union[sp_sparse.csr_matrix, np.ndarray]],
        batch_indices_per_label: List[Union[List[int], np.ndarray]] = None,
        gene_names: Union[List[str], np.ndarray] = None,
    ):
        """Populates the data attributes of a GeneExpressionDataset object from a ``n_labels``-long
            ``list`` of (nb_cells, nb_genes) matrices.

        :param Xs: RNA counts in the form of a list of np.ndarray with shape (..., nb_genes)
        :param batch_indices_per_label: cell-wise batch indices, for each cell label.
        :param gene_names: gene names, stored as ``str``.
        """
        nb_genes = Xs[0].shape[1]
        if not all(X.shape[1] == nb_genes for X in Xs):
            raise ValueError("All batches must have same nb_genes")

        X = np.concatenate(Xs) if type(Xs[0]) is np.ndarray else sp_sparse.vstack(Xs)
        labels = np.concatenate(
            [i * np.ones(cluster.shape[0], dtype=np.int64) for i, cluster in enumerate(Xs)]
        )
        batch_indices = (
            np.concatenate(batch_indices_per_label).astype(np.int64)
            if batch_indices_per_label
            else None
        )

        self.populate_from_data(
            X=X, batch_indices=batch_indices, labels=labels, gene_names=gene_names
        )

    def populate_from_datasets(
        self,
        gene_datasets_list: List["GeneExpressionDataset"],
        on: str = "gene_names",
        sharing_intstructions_dict: Dict[str, bool] = None,
    ):
        """Populates the data attribute of a GeneExpressionDataset
        from multiple ``GeneExperessionDataset`` objects, merged using the intersection
        of a gene-wise attribute (``gene_names`` by default).
        Note that datasets should all have gene_dataset.n_labels=0.
        Batch indices are generated in the same order as datasets are given.

        :param gene_datasets: a sequence of ``GeneExpressionDataset`` objects.
        :param on: attribute to select gene interesection
        :param sharing_intstructions_dict: Instructions on how to share cell-wise attributes between datasets.
            Keys are the attribute name and values are either:
                * "offset": to add an offset corresponding to the number of categories already existing.
                e.g for batch_indices, if the first dataset has batches 0 and 1,
                in the merged dataset, the second dataset's batch indices will start at 2.
                * "concatenate": concatenate the attibute, no changes applied.
        """
        # sanity check
        if not all([hasattr(gene_dataset, on) for gene_dataset in gene_datasets_list]):
            raise ValueError(
                "All datasets should have the merge key 'on' as an attribute"
            )

        # set default sharing behaviour
        if sharing_intstructions_dict is None:
            sharing_intstructions_dict = {"batch_indices": "offset"}

        # get insterection based on gene attribute `on` and get attribute intersection
        genes_to_keep = set.intersection(
            *[set(getattr(gene_dataset, on)) for gene_dataset in gene_datasets_list]
        )
        gene_attributes_to_keep = set.intersection(
            *[
                set(gene_dataset.gene_attribute_names)
                for gene_dataset in gene_datasets_list
            ]
        )
        cell_attributes_to_keep = set.intersection(
            *[
                set(gene_dataset.cell_attribute_names)
                for gene_dataset in gene_datasets_list
            ]
        )

        # keep gene order and attributes of the first dataset
        gene_to_keep = [
            gene for gene in getattr(gene_datasets_list[0], on) if gene in genes_to_keep
        ]
        logger.info("Keeping {nb_genes} genes".format(nb_genes=len(genes_to_keep)))

        # filter genes
        Xs = []
        for dataset in gene_datasets_list:
            dataset.filter_genes(gene_to_keep, on=on)
            Xs.append(dataset.X)

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

        self.populate_from_data(X=X)

        # keep gene attributes of first dataset, and keep all mappings (e.g gene types)
        for attribute_name in gene_attributes_to_keep:
            self.initialize_gene_attribute(
                attribute_name, getattr(gene_datasets_list[0], attribute_name)
            )
            for gene_dataset in gene_datasets_list:
                mapping_names = gene_dataset.attribute_mappings[attribute_name]
                for mapping_name in mapping_names:
                    self.initialize_mapped_attribute(
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
                        for gene_dataset in gene_datasets_list
                    ]
                )
            )

            attribute_values = []

            if instruction == "concatenate":
                mappings = defaultdict(set)
                # create new mapping
                for mapping_name in mapping_names_to_keep:
                    mappings[mapping_name] = mappings[mapping_name].union(
                        getattr(self, mapping_name)
                    )
                mappings = {k: list(v) for k, v in mappings.items()}

                for gene_dataset in gene_datasets_list:
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
                for i, gene_dataset in enumerate(gene_datasets_list):
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

            self.initialize_cell_attribute(
                attribute_name, np.concatenate(attribute_values)
            )

            for mapping_name, mapping_values in mappings.items():
                self.initialize_mapped_attribute(
                    attribute_name, mapping_name, mapping_values
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
        logger.info("Computing the library size for the new data")
        self.compute_library_size_batch()

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
        """Sets batch indices and computes the n_batches"""
        self.n_batches = len(np.unique(batch_indices))
        self._batch_indices = batch_indices

    @property
    def labels(self) -> np.ndarray:
        return self._labels

    @labels.setter
    def labels(self, labels: Union[List[int], np.ndarray]):
        """Sets labels and computes n_labels"""
        self.n_labels = len(np.unique(labels))
        self._labels = labels

    def remap_cell_types(self, labels):
        """Remaps cell_types using new labels."""
        new_cell_types = []
        n_unknown_cell_types = 0
        for new_label in np.unique(labels).astype(np.uint16):
            if new_label < self.n_labels:
                new_cell_types.append(getattr(self, "cell_types")[new_label])
            # if new cell_type, needs to be set elsewhere, using 'unknown_c...' in the meantime
            else:
                new_cell_types.append("unknown_cell_type_" + str(n_unknown_cell_types))
                n_unknown_cell_types += 1

    @property
    def norm_X(self) -> Union[sp_sparse.csr_matrix, np.ndarray]:
        """Returns a normalized version of X."""
        return self._norm_X

    @norm_X.setter
    def norm_X(self, norm_X: Union[sp_sparse.csr_matrix, np.ndarray]):
        self._norm_X = norm_X
        self.register_dataset_version("norm_X")

    @property
    def corrupted_X(self) -> Union[sp_sparse.csr_matrix, np.ndarray]:
        """Returns the corrupted version of X."""
        return self._corrupted_X

    @corrupted_X.setter
    def corrupted_X(self, corrupted_X: Union[sp_sparse.csr_matrix, np.ndarray]):
        self._corrupted_X = corrupted_X
        self.register_dataset_version("corrupted_X")

    def remap_categorical_attributes(self):
        for attribute_name in self.cell_categorical_attribute_names:
            logger.info("Remapping %s to [0,N]" % attribute_name)
            attr = getattr(self, attribute_name)
            new_attr, _ = remap_categories(attr)
            setattr(self, attribute_name, new_attr)

    def register_dataset_version(self, version_name):
        """Registers a version of the dataset, e.g normalized version."""
        self.dataset_versions.add(version_name)

    def initialize_cell_attribute(self, attribute_name, attribute, categorical=False):
        """Sets and registers a cell-wise attribute, e.g annotation information."""
        if not self.nb_cells == len(attribute):
            raise ValueError(
                "Number of cells ({n_cells}) and length of cell attribute ({n_attr}) mismatch".format(
                    n_cells=self.nb_cells, n_attr=len(attribute)
                )
            )
        setattr(self, attribute_name, attribute)
        self.cell_attribute_names.add(attribute_name)
        if categorical:
            self.cell_categorical_attribute_names.add(attribute_name)

    def initialize_gene_attribute(self, attribute_name, attribute):
        """Sets and registers a gene-wise attribute, e.g annotation information."""
        if not self.nb_genes == len(attribute):
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
        setattr(self, mapping_name, mapping_values)

    def compute_library_size_batch(self):
        """Computes the library size per batch."""
        self.local_means = np.zeros((self.nb_cells, self.nb_genes))
        self.local_vars = np.zeros((self.nb_cells, self.nb_genes))
        for i_batch in range(self.n_batches):
            idx_batch = (self.batch_indices == i_batch).ravel()
            self.local_means[idx_batch], self.local_vars[
                idx_batch
            ] = compute_library_size(self.X[idx_batch])
        self.cell_attribute_names.update(["local_means", "local_vars"])

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
        return (
            X_batch,
            torch.from_numpy(self.local_means[indices].astype(np.float32)),
            torch.from_numpy(self.local_vars[indices].astype(np.float32)),
            torch.from_numpy(self.batch_indices[indices].astype(np.int64)),
            torch.from_numpy(self.labels[indices].astype(np.int64)),
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
            if subset_genes.dtype is not np.bool
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
        logger.info("Filtering non-expressing cells.")
        self.filter_cells(1)

    def filter_cells(self, min_count=1):
        mask_cells_to_keep = np.asarray(self.X.sum(axis=1) >= min_count)
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

        # update cell-related attributes accordingly
        for attribute_name in self.cell_attribute_names:
            attr = getattr(self, attribute_name)
            setattr(self, attribute_name, attr[subset_cells])

        # subsample datasets which calls library_size_batch which requires already updated attributes
        self.X = self.X[subset_cells]
        for version_name in self.dataset_versions:
            data = getattr(self, version_name)
            setattr(self, version_name, data[:, subset_cells])

        logger.info(
            "Downsampled from {nb_cells} to {new_nb_cells} cells".format(
                nb_cells=nb_cells_old, new_nb_cells=self.nb_cells
            )
        )

    def subsample_genes(
        self,
        new_n_genes: Optional[int] = None,
        new_ratio_genes: Optional[float] = None,
        subset_genes: Optional[Union[List[int], List[bool], np.ndarray]] = None,
    ):
        """Wrapper around ``update_genes`` allowing for manual and automatic (based on count variance) subsampling.

        The function either:
            * Subsamples `new_n_genes` genes among all genes
            * Subsambles a proportion of `new_ratio_genes` of the genes
            * Subsamples the genes in `subset_genes`

        :param subset_genes: list of indices or mask of genes to retain
        :param new_n_genes: number of genes to retain, the highly variable genes will be kept
        :param new_ratio_genes: proportion of genes to retain, the highly variable genes will be kept
        """
        if new_ratio_genes is not None:
            if 0 < new_ratio_genes < 1:
                new_n_genes = int(new_ratio_genes * self.nb_genes)
            else:
                logger.info(
                    "Not subsampling. Expecting 0 < (new_ratio_genes={new_ratio_genes})  < 1.".format(
                        new_ratio_genes=new_ratio_genes
                    )
                )
                return

        if new_n_genes is not None:
            if new_n_genes >= self.nb_genes or new_n_genes < 1:
                logger.info(
                    "Not subsampling. Expecting: 1 < (new_n_genes={new_n_genes}) <= self.nb_genes".format(
                        new_n_genes=new_n_genes
                    )
                )
                return

            std_scaler = StandardScaler(with_mean=False)
            std_scaler.fit(self.X.astype(np.float64))
            subset_genes = np.argsort(std_scaler.var_)[::-1][:new_n_genes]

        if subset_genes is None:
            logger.info(
                "Not subsampling. No parameter given".format(new_n_genes=new_n_genes)
            )
            return

        self.update_genes(np.array(subset_genes))

    def subsample_cells(self, size: Union[int, float] = 1.0):
        """Wrapper around ``update_cells`` allowing for automatic (based on sum of counts) subsampling.
        If size is a:
            * (0,1) float: subsample 100*``size`` % of the cells
            * int: subsample ``size`` cells
        """
        new_n_cells = int(size * self.nb_genes) if type(size) is not int else size
        indices = np.argsort(np.asarray(self.X.sum(axis=1)).ravel())[::-1][:new_n_cells]
        self.update_cells(indices)

    def reorder_genes(self, first_genes: Union[List[str], np.ndarray]):
        """
        Performs a in-place reordering of genes and gene-related attributes.

        Reorder genes according to the ``first_genes`` list of gene names.
        Consequently, modifies in-place the data ``X`` and the registered gene attributes.

        :param first_genes: New ordering of the genes; if some genes are missing, they will be added after the
                            first_genes in the same order as they before
        """

        _, _, new_order_first = np.intersect1d(
            first_genes, self.gene_names, return_indices=True
        )
        new_order_second = [x for x in range(self.nb_genes) if x not in new_order_first]
        new_order = np.hstack([new_order_first, new_order_second])

        # update datasets
        self.X = self.X[:, new_order]
        for version_name in self.dataset_versions:
            data = getattr(self, version_name)
            setattr(self, version_name, data[:, new_order])

        # update gene-related attributes accordingly
        for attribute_name in self.gene_attribute_names:
            attr = getattr(self, attribute_name)
            setattr(self, attribute_name, attr[new_order])

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
            np.where(getattr(self, "cell_types") == cell_type)[0][0]
            for cell_type in cell_types
        ]
        return np.asarray(labels, dtype=np.int64)

    def _gene_idx(self, genes):
        if type(genes[0]) is not int:
            genes_idx = [
                np.where(gene == getattr(self, "gene_names"))[0][0] for gene in genes
            ]
        else:
            genes_idx = genes
        return np.asarray(genes_idx, dtype=np.int64)

    def filter_cell_types(self, cell_types: Union[List[str], List[int], np.ndarray]):
        """Performs in-place filtering of cells by keeping cell types in ``cell_types``.

        :param cell_types: numpy array of type np.int (indices) or np.str (cell-types names)
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
        cell_types: Union[
            Tuple[int, ...], Tuple[str, ...], List[int], List[str], np.ndarray
        ],
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

    def normalize(self):
        scaling_factor = self.X.mean(axis=1)
        self.norm_X = self.X / scaling_factor.reshape(len(scaling_factor), 1)

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
        self.corrupted_X = copy.deepcopy(self.X)
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
        ):  # replace the entry n with a Bin(n, 0.2) random variable.
            i, j = (k.ravel() for k in np.indices(self.X.shape))
            ix = np.random.choice(len(i), int(np.floor(rate * len(i))), replace=False)
            i, j = i[ix], j[ix]
            self.corrupted_X[i, j] = np.random.binomial(
                n=(self.X[i, j]).astype(np.int32), p=0.2
            )
        else:
            raise NotImplementedError("Unknown corruption method.")

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

    # check lengths, allow absent cell types
    if n_categories > len(mapping_from):
        raise ValueError(
            "Size of provided mapping_from greater than the number of categories."
        )
    if len(mapping_to) != len(mapping_from):
        raise ValueError("Length mismatch between mapping_to and mapping_from.")

    new_categories = np.copy(original_categories)
    for idx_from, idx_to in zip(mapping_from, mapping_to):
        new_categories[original_categories == idx_from] = idx_to
    return new_categories.astype(np.uint16), n_categories


def compute_library_size(
    data: Union[sp_sparse.csr_matrix, np.ndarray]
) -> Tuple[np.ndarray, np.ndarray]:
    log_counts = np.log(data.sum(axis=1))
    local_mean = (np.mean(log_counts) * np.ones((data.shape[0], 1))).astype(np.float32)
    local_var = (np.var(log_counts) * np.ones((data.shape[0], 1))).astype(np.float32)
    return local_mean, local_var


class SpatialGeneExpressionDataset(GeneExpressionDataset):
    """Generic class representing RNA counts with spatial coordinates.


    :param X: RNA counts matrix, sparse format supported (e.g ``scipy.sparse.csr_matrix``).
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
        gene_names: Union[List[str], np.ndarray, sp_sparse.csr_matrix] = None,
        cell_types: Union[List[int], np.ndarray] = None,
        x_coord: Union[List[float], np.ndarray, sp_sparse.csr_matrix] = None,
        y_coord: Union[List[float], np.ndarray, sp_sparse.csr_matrix] = None,
    ):
        super().__init__()
        super().populate_from_data(
            X=X,
            batch_indices=batch_indices,
            labels=labels,
            gene_names=gene_names,
            cell_types=cell_types,
            cell_attributes_dict={"x_coord": x_coord, "y_coord": y_coord},
        )

    def make_tensor_batch_from_indices(
        self, X_batch: Union[sp_sparse.csr_matrix, np.ndarray], indices: np.ndarray
    ) -> Tuple[torch.Tensor, ...]:
        """Given indices and batch X_batch, returns a full batch of ``Torch.Tensor``"""
        if isinstance(X_batch, np.ndarray):
            X_batch = torch.from_numpy(X_batch)
        else:
            X_batch = torch.from_numpy(X_batch.toarray().astype(np.float32))
        return (
            X_batch,
            torch.from_numpy(self.local_means[indices].astype(np.float32)),
            torch.from_numpy(self.local_vars[indices].astype(np.float32)),
            torch.from_numpy(self.batch_indices[indices].astype(np.int64)),
            torch.from_numpy(self.labels[indices].astype(np.int64)),
            torch.from_numpy(getattr(self, "x_coord")[indices].astype(np.float32)),
            torch.from_numpy(getattr(self, "y_coord")[indices].astype(np.float32)),
        )


class DownloadableDataset(GeneExpressionDataset, ABC):
    """Sub-class of ``GeneExpressionDataset`` which downloads its data to disk and
    then populates its attributes with it.

    In particular, it has a ``delayed_populating`` parameter allowing for instantiation
    without populating the attributes.

    :param urls: single or multiple url.s from which to download the data.
    :param filenames: filenames for the downloaded data.
    :param save_path: path to data storage.
    :param delayed_populating: If False, populate object upon isntantiation.
        Else, allow for a delayed manual call to ``populate`` method.
    """

    def __init__(
        self,
        urls: Union[str, List[str]],
        filenames: Union[str, List[str]] = None,
        save_path: str = "data/",
        delayed_populating: bool = False,
    ):
        super().__init__()
        if isinstance(urls, str):
            self.urls = [urls]
        elif urls is None:
            self.urls = []
        else:
            self.urls = urls
        if isinstance(filenames, str):
            self.filenames = [filenames]
        elif filenames is None:
            self.filenames = ["dataset_{i}".format(i=i) for i in range(len(self.urls))]
        else:
            self.filenames = filenames

        self.save_path = save_path
        self.download()
        if not delayed_populating:
            self.populate()

    def download(self):
        for url, download_name in zip(self.urls, self.filenames):
            _download(url, self.save_path, download_name)

    @abstractmethod
    def populate(self):
        """Populates a ``DonwloadableDataset`` object's data attributes.

        E.g by calling one of ``GeneExpressionDataset``'s ``populate_from...`` methods.
        """
        pass


def _download(url, save_path, filename):
    """Writes data from url to file."""
    if os.path.exists(os.path.join(save_path, filename)):
        logger.info("File %s already downloaded" % (os.path.join(save_path, filename)))
        return

    r = urllib.request.urlopen(url)
    logger.info("Downloading file at %s" % os.path.join(save_path, filename))

    def read_iter(file, block_size=1000):
        """Given a file 'file', returns an iterator that returns bytes of
        size 'blocksize' from the file, using read()."""
        while True:
            block = file.read(block_size)
            if not block:
                break
            yield block

    # Create the path to save the data
    if not os.path.exists(save_path):
        os.makedirs(save_path)

    with open(os.path.join(save_path, filename), "wb") as f:
        for data in read_iter(r):
            f.write(data)
