import csv
import logging
import os

from typing import Dict, List, Union

import numpy as np
import scipy.sparse as sp_sparse

from scvi.dataset import DownloadableDataset

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
        save_path="data/",
        genes_to_keep=None,
        genes_fish=None,
        additional_genes=558,
        delayed_populating: bool = False,
    ):
        # handle mutable defaults
        genes_to_keep = genes_to_keep if genes_to_keep else []
        genes_fish = genes_fish if genes_fish else []

        # If we want to harmonize the dataset with the OsmFISH dataset, we need to keep
        # OsmFISH genes and order the genes from Cortex accordingly
        self.genes_fish = genes_fish
        # If there are specific genes we'd like to keep
        self.genes_to_keep = genes_to_keep
        # Number of genes we want to keep
        self.additional_genes = additional_genes

        super().__init__(
            urls="https://storage.googleapis.com/linnarsson-lab-www-blobs/blobs"
            "/cortex/expression_mRNA_17-Aug-2014.txt",
            filenames="expression.bin",
            save_path=save_path,
            delayed_populating=delayed_populating,
        )

    def load_from_disk(self):
        logger.info("Loading Cortex data")
        rows = []
        gene_names = []
        with open(os.path.join(self.save_path, self.filenames[0]), "r") as csvfile:
            data_reader = csv.reader(csvfile, delimiter="\t")
            clusters = None
            for i, row in enumerate(data_reader):
                if i == 1:
                    precise_clusters = np.array(row, dtype=str)[2:]
                if i == 8:  # 7 + 1 in pandas
                    clusters = np.asarray(row, dtype=str)[2:]
                if i >= 11:  # 10 + 1 in pandas
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
        additional_genes = []
        for gene_cortex in range(len(gene_names)):
            for gene_fish in self.genes_fish:
                if gene_names[gene_cortex].lower() == gene_fish.lower():
                    additional_genes.append(gene_cortex)
        for gene_cortex in range(len(gene_names)):
            for gene_fish in self.genes_to_keep:
                if gene_names[gene_cortex].lower() == gene_fish.lower():
                    additional_genes.append(gene_cortex)

        selected = np.std(X, axis=0).argsort()[-self.additional_genes:][::-1]
        selected = np.unique(np.concatenate((selected, np.array(additional_genes))))
        selected = np.array([int(select) for select in selected])
        X = X[:, selected]
        gene_names = gene_names[selected]

        # Then we reorganize the genes so that the genes from the smFISH dataset
        # appear first
        if len(self.genes_fish) > 0:
            X, gene_names = reorder_genes(X, gene_names, self.genes_fish)
            umi = np.sum(X[:, :len(self.genes_fish)], axis=1)
            X = X[umi > 10, :]
            labels = labels[umi > 10]

        logging.info("Finished preprocessing Cortex data")
        return {
            "X": X,
            "labels": labels,
            "gene_names": gene_names,
            "cell_types": cell_types,
            "cell_attributes_dict": cell_attributes_dict
        }

    def instantiate_gene_expression_dataset(
        self,
        X: Union[np.ndarray, sp_sparse.csr_matrix],
        labels: Union[List[int], np.ndarray, sp_sparse.csr_matrix] = None,
        gene_names: Union[List[str], np.ndarray, sp_sparse.csr_matrix] = None,
        cell_types: Union[List[int], np.ndarray] = None,
        cell_attributes_dict: Dict[str, Union[List, np.ndarray]] = None,
    ):
        super(DownloadableDataset, self).__init__(
            X=X,
            labels=labels,
            gene_names=gene_names,
            cell_types=cell_types,
            cell_attributes_dict=cell_attributes_dict,
        )


