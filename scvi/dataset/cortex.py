import csv
import logging
import os

from typing import Dict, List, Union, Optional

import numpy as np
import scipy.sparse as sp_sparse

from scvi.dataset.dataset import DownloadableDataset

logger = logging.getLogger(__name__)


class CortexDataset(DownloadableDataset):
    r""" Loads cortex dataset.

    The `Mouse Cortex Cells dataset`_ contains 3005 mouse cortex cells and gold-standard labels for
    seven distinct cell types. Each cell type corresponds to a cluster to recover. We retain top 558 genes
    ordered by variance.

    Args:
        :save_path: Save path of raw data file. Default: ``'data/'``.

    Examples:
        >>> gene_dataset = CortexDataset()

    .. _Mouse Cortex Cells dataset:
        https://storage.googleapis.com/linnarsson-lab-www-blobs/blobs/cortex/expression_mRNA_17-Aug-2014.txt

    """

    def __init__(
        self,
        save_path: str = "data/",
        genes_to_keep: Optional[List[str]] = None,
        total_genes: Optional[int] = None,
        delayed_populating: bool = False,
    ):
        self.genes_to_keep = genes_to_keep
        self.total_genes = total_genes

        self.precise_labels = None

        super().__init__(
            urls="https://storage.googleapis.com/linnarsson-lab-www-blobs/blobs"
                 "/cortex/expression_mRNA_17-Aug-2014.txt",
            filenames="expression.bin",
            save_path=save_path,
            delayed_populating=delayed_populating,
        )

    def populate(self):
        data = self.load_from_disk()
        data = self.preprocess(**data)
        self.populate_from_data(**data)

    def load_from_disk(self):
        logger.info("Loading Cortex data")
        rows = []
        gene_names = []
        with open(os.path.join(self.save_path, self.filenames[0]), "r") as csvfile:
            data_reader = csv.reader(csvfile, delimiter="\t")
            for i, row in enumerate(data_reader):
                if i == 1:
                    precise_clusters = np.array(row, dtype=str)[2:]
                if i == 8:
                    clusters = np.asarray(row, dtype=str)[2:]
                if i >= 11:
                    rows.append(row[1:])
                    gene_names.append(row[0])
        cell_types, labels = np.unique(clusters, return_inverse=True)
        _, self.precise_labels = np.unique(precise_clusters, return_inverse=True)
        X = np.array(rows, dtype=np.int).T[1:]
        gene_names = np.array(gene_names, dtype=np.str)
        return {
            "X": X,
            "labels": labels,
            "gene_names": gene_names,
            "cell_types": cell_types,
            "cell_attributes_dict": {"precise_labels": precise_clusters},
        }

    def preprocess(
        self,
        X: Union[np.ndarray, sp_sparse.csr_matrix],
        labels: Union[List[int], np.ndarray, sp_sparse.csr_matrix] = None,
        gene_names: Union[List[str], np.ndarray, sp_sparse.csr_matrix] = None,
        cell_types: Union[List[int], np.ndarray] = None,
        cell_attributes_dict: Dict[str, Union[List, np.ndarray]] = None,
    ):
        gene_indices = []
        if self.genes_to_keep is not None:
            _, gene_indices, _ = np.intersect1d(
                self.genes_to_keep, gene_names, return_indices=True
            )
            gene_indices = list(gene_indices)

        if self.total_genes is not None:
            genes_by_variance = np.std(X, axis=0).argsort()[::-1]
            to_add = self.total_genes - len(gene_indices)
            added = set(gene_indices)

            for gene_index in genes_by_variance:
                if to_add == 0:
                    break
                if gene_index not in added:
                    added.add(gene_index)
                    gene_indices.append(gene_index)
                    to_add -= 1

        gene_indices = np.array(gene_indices)
        if self.total_genes is None and self.genes_to_keep is None:
            gene_indices = slice(None)

        X = X[:, gene_indices]
        gene_names = gene_names[gene_indices]

        logging.info("Finished preprocessing Cortex data")
        return {
            "X": X,
            "labels": labels,
            "gene_names": gene_names,
            "cell_types": cell_types,
            "cell_attributes_dict": cell_attributes_dict
        }

