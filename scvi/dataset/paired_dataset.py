# -*- coding: utf-8 -*-

"""Handling datasets.
For the moment, is initialized with a torch Tensor of size (n_cells, nb_genes)"""
import logging

from typing import Dict, List, Tuple, Union, Optional

import numpy as np
import scipy.sparse as sp_sparse
import torch
from scvi.dataset.dataset import GeneExpressionDataset

logger = logging.getLogger(__name__)


class PairedExpressionDataset(GeneExpressionDataset):
    """Generic class representing RNA counts and a paired measurement.

    This class inherits scVI's base dataset class. See docstring of
    GeneExpressionDataset for more information.
    """

    def __init__(self):
        super().__init__()
        self._Y = None
        self.paired_names = None
        self.paired_attribute_names = set()

    def populate_paired_measurements(self, Y, paired_names: Union[List[str], np.ndarray] = None):
        self._Y = (
            np.ascontiguousarray(Y, dtype=np.float32)
            if isinstance(Y, np.ndarray)
            else Y
        )
        if paired_names is not None:
            self.initialize_paired_attribute(
                "paired_names", np.char.lower(np.asarray(paired_names, dtype=np.str))
            )

    def populate_from_data(
        self,
        X: Union[np.ndarray, sp_sparse.csr_matrix],
        Y: Union[np.ndarray, sp_sparse.csr_matrix],
        batch_indices: Union[List[int], np.ndarray, sp_sparse.csr_matrix] = None,
        labels: Union[List[int], np.ndarray, sp_sparse.csr_matrix] = None,
        gene_names: Union[List[str], np.ndarray] = None,
        paired_names: Union[List[str], np.ndarray] = None,
        cell_types: Union[List[int], np.ndarray] = None,
        cell_attributes_dict: Dict[str, Union[List, np.ndarray]] = None,
        gene_attributes_dict: Dict[str, Union[List, np.ndarray]] = None,
        paired_attributes_dict: Dict[str, Union[List, np.ndarray]] = None,
    ):
        """Populates the data attributes of a PairedExpressionDataset object.
        
        This takes two matirces. One is a (nb_cells, nb_genes) matrix and the other is a (nb_cells, nb_paired) matrix.
        This assumes cells are from the two matrices are already in the same order.

        :param X: RNA counts matrix, sparse format supported (e.g ``scipy.sparse.csr_matrix``).
        :param Y: Paired count matrix (e.g CITE-seq protein measurements), sparse format supported
            (e.g ``scipy.sparse.csr_matrix``).
        :param batch_indices: ``np.ndarray`` with shape (nb_cells,). Maps each cell to the batch
            it originates from. Note that a batch most likely refers to a specific piece
            of tissue or a specific experimental protocol.
        :param labels: ``np.ndarray`` with shape (nb_cells,). Cell-wise labels. Can be mapped
            to cell types using attribute mappings.
        :param gene_names: ``List`` or ``np.ndarray`` with length/shape (nb_genes,).
            Maps each gene to its name.
        :param paired_names: ``List`` or ``np.ndarray`` with length/shape (Y.shape[1],).
            Maps each paired measurement to its name.
        :param cell_types: Maps each integer label in ``labels`` to a cell type.
        :param cell_attributes_dict: ``List`` or ``np.ndarray`` with shape (nb_cells,).
        :param gene_attributes_dict: ``List`` or ``np.ndarray`` with shape (nb_genes,).
        :param paired_attributes_dict: ``List`` or ``np.ndarray`` with shape (nb_paired,).
        """

        super().populate_from_data(
            X,
            batch_indices,
            labels,
            gene_names,
            cell_types,
            cell_attributes_dict,
            gene_attributes_dict,
        )

        self.populate_paired_measurements(Y, paired_names)

    def populate_from_per_batch_array(
        self,
        X: np.ndarray,
        Y: np.ndarray,
        labels_per_batch: Union[np.ndarray, List[Union[List[int], np.ndarray]]] = None,
        gene_names: Union[List[str], np.ndarray] = None,
        paired_names: Union[List[str], np.ndarray] = None,
    ):
        """Populates the data attributes of a PairedExpressionDataset object
        from an array with shape (n_batches, nb_cells, nb_genes) and an array
        with shape (n_batches, nb_cells, nb_paired)

        :param X: np.ndarray with shape (n_batches, nb_cells, nb_genes).
        :param Y: np.ndarray with shape (n_batches, nb_cells, nb_paired).
        :param labels_per_batch: cell-wise labels for each batch.
        :param gene_names: gene names, stored as ``str``.
        :param paired_names: paired names, stored as ``str``.
        """

        super().populate_from_per_batch_array(
            X,
            labels_per_batch,
            gene_names,
        )
        if len(np.squeeze(X.shape)) != 3:
            raise ValueError(
                "Shape of np.squeeze(Y) != 3. Use populate_from_data "
                "if your dataset has shape (nb_cells, nb_paired)"
            )
        Y = Y.reshape(-1, Y.shape[2])
        self.populate_paired_measurements(Y, paired_names)

    def populate_from_per_batch_list(
        self,
        Xs: np.ndarray,
        Ys: np.ndarray,
        labels_per_batch: Union[np.ndarray, List[Union[List[int], np.ndarray]]] = None,
        gene_names: Union[List[str], np.ndarray] = None,
        paired_names: Union[List[str], np.ndarray] = None,
    ):
        """Populates the data attributes of a PairedExpressionDataset object from a ``n_batches``-long
            ``list`` of (nb_cells, nb_genes) matrices paired with ``list`` of (nb_cells, nb_paired) matrices

        :param Xs: RNA counts in the form of a list of np.ndarray with shape (..., nb_genes)
        :param Ys: Paired counts in the form of a list of np.ndarray with shape (..., nb_paired)
        :param labels_per_batch: list of cell-wise labels for each batch.
        :param gene_names: gene names, stored as ``str``.
        :param paired_names: paired names, stored as ``str``.
        """
        super().populate_from_per_batch_list(
            Xs,
            labels_per_batch,
            gene_names,
        )
        nb_paired = Ys[0].shape[1]
        if not all(Y.shape[1] == nb_paired for Y in Ys):
            raise ValueError("All batches must have same nb_paired")
        Y = np.concatenate(Ys) if type(Ys[0]) is np.ndarray else sp_sparse.vstack(Ys)
        self.populate_paired_measurements(Y, paired_names)

    def populate_from_per_label_list(
        self,
        Xs: List[Union[sp_sparse.csr_matrix, np.ndarray]],
        Ys: List[Union[sp_sparse.csr_matrix, np.ndarray]],
        batch_indices_per_label: List[Union[List[int], np.ndarray]] = None,
        gene_names: Union[List[str], np.ndarray] = None,
        paired_names: Union[List[str], np.ndarray] = None,
    ):
        """Populates the data attributes of a PairedExpressionDataset object from a ``n_labels``-long
            ``list`` of (nb_cells, nb_genes) matrices paired with a ``list`` of (nb_cells, nb_paired) matrices.

        :param Xs: RNA counts in the form of a list of np.ndarray with shape (..., nb_genes)
        :param Ys: Paired counts in the form of a list of np.ndarray with shape (..., nb_paired)
        :param batch_indices_per_label: cell-wise batch indices, for each cell label.
        :param gene_names: gene names, stored as ``str``.
        :param paired_names: paired names, stored as ``str``.
        """
        super().populate_from_per_label_list(Xs, batch_indices_per_label. gene_names)
        nb_paired = Ys[0].shape[1]
        if not all(Y.shape[1] == nb_paired for Y in Ys):
            raise ValueError("All batches must have same nb_paired")

        Y = np.concatenate(Ys) if type(Ys[0]) is np.ndarray else sp_sparse.vstack(Ys)
        self.populate_paired_measurements(Y, paired_names)

    def populate_from_datasets(
        self,
        gene_datasets_list: List["PairedExpressionDataset"],
        on: str = "gene_names",
        sharing_intstructions_dict: Dict[str, bool] = None,
    ):
        raise NotImplementedError

    @property
    def Y(self):
        return self._Y

    @Y.setter
    def Y(self, Y: Union[np.ndarray, sp_sparse.csr_matrix]):
        """Sets the data attribute ``Y``."""
        self._Y = Y

    @property
    def nb_paired(self) -> int:
        return self.Y.shape[1]

    def initialize_paired_attribute(self, attribute_name, attribute):
        """Sets and registers a column-wise attribute, e.g protein name."""
        if not self.nb_paired == len(attribute):
            raise ValueError(
                "Number of paired measurements ({n_paired}) and length of paired attribute ({n_attr}) mismatch".format(
                    n_paired=self.nb_paired, n_attr=len(attribute)
                )
            )
        setattr(self, attribute_name, attribute)
        self.paired_attribute_names.add(attribute_name)

    def collate_fn(
        self, batch: Union[List[int], np.ndarray]
    ) -> Tuple[torch.Tensor, ...]:
        """Batch creation function to be passed to Torch's DataLoader."""
        indices = np.asarray(batch)
        X = self.X[indices]
        Y = self.Y[indices]
        return self.make_tensor_batch_from_indices(X, Y, indices)

    def collate_fn_corrupted(
        self, batch: Union[List[int], np.ndarray]
    ) -> Tuple[torch.Tensor, ...]:
        """Batch creation function to be passed to Torch's DataLoader."""
        indices = np.asarray(batch)
        X_batch = self.corrupted_X[indices]
        Y_batch = self.Y[indices]
        return self.make_tensor_batch_from_indices(X_batch, Y_batch, indices)

    def make_tensor_batch_from_indices(
        self, X_batch: Union[sp_sparse.csr_matrix, np.ndarray], Y_batch: Union[sp_sparse.csr_matrix, np.ndarray], indices: np.ndarray
    ) -> Tuple[torch.Tensor, ...]:
        """Given indices and batches X_batch, Y_batch, returns a full batch of ``Torch.Tensor``"""
        if isinstance(X_batch, np.ndarray):
            X_batch = torch.from_numpy(X_batch)
        else:
            X_batch = torch.from_numpy(X_batch.toarray().astype(np.float32))
        if isinstance(Y_batch, np.ndarray):
            Y_batch = torch.from_numpy(Y_batch)
        else:
            Y_batch = torch.from_numpy(Y_batch.toarray().astype(np.float32))
        return (
            X_batch,
            Y_batch,
            torch.from_numpy(self.local_means[indices].astype(np.float32)),
            torch.from_numpy(self.local_vars[indices].astype(np.float32)),
            torch.from_numpy(self.batch_indices[indices].astype(np.int64)),
            torch.from_numpy(self.labels[indices].astype(np.int64)),
        )

    def update_cells(self, subset_cells):
        """Performs a in-place sub-sampling of cells and cell-related attributes.

        Sub-selects cells according to ``subset_cells`` sub-index.
        Consequently, modifies in-place the data ``X`` and `Y``, its versions and the registered cell attributes.

        :param subset_cells: Index used for cell sub-sampling.
            Either a ``int`` array with arbitrary shape which values are the indexes of the cells to keep.
            Or boolean array used as a mask-like index.
        """
        self.Y = self.Y[subset_cells]
        super().update_cells(subset_cells)

