from unittest import TestCase

from scvi.dataset import CortexDataset
from .utils import unsupervised_training_one_epoch
import numpy as np


class TestCortexDataset(TestCase):
    def test_populate(self):
        dataset = CortexDataset(save_path="tests/data")
        unsupervised_training_one_epoch(dataset)

    def test_variance_and_order_and_size(self):
        to_keep = ["THY1", "SST", "TOMEM2", "MOG"]
        total_genes = 10
        dataset_full = CortexDataset(save_path="tests/data", total_genes=None)
        dataset_small = CortexDataset(
            save_path="tests/data", genes_to_keep=to_keep, total_genes=total_genes
        )
        self.assertListEqual(dataset_small.gene_names[:4].tolist(), to_keep)

        small_variance = np.std(dataset_small.X[:, 4:], axis=0).argsort()[::-1]
        self.assertListEqual(small_variance.tolist(), list(range(6)))

        full_variance = np.std(dataset_full.X, axis=0).argsort()[::-1]
        variable_genes_all = dataset_full.gene_names[full_variance]
        genes_truth = (to_keep + [g for g in variable_genes_all if g not in to_keep])[
            :total_genes
        ]
        self.assertListEqual(dataset_small.gene_names.tolist(), genes_truth)
