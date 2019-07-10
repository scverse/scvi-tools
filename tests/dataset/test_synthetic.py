from unittest import TestCase

import numpy as np

from scvi.dataset import (
    SyntheticDataset,
    SyntheticRandomDataset,
    SyntheticDatasetCorr,
    ZISyntheticDatasetCorr,
)
from .utils import unsupervised_training_one_epoch


class TestSyntheticDataset(TestCase):
    def test_train_one(self):
        dataset = SyntheticDataset(batch_size=10, nb_genes=10)
        unsupervised_training_one_epoch(dataset)

    def test_RandomDataset_populate_and_train_one(self):
        dataset = SyntheticRandomDataset(save_path="tests/data")
        unsupervised_training_one_epoch(dataset)

    def test_DatasetCorr_populate_and_train_one(self):
        dataset = SyntheticDatasetCorr(n_cells_cluster=10)
        self.assertListEqual(
            np.unique(dataset.labels).tolist(), np.arange(dataset.n_clusters).tolist()
        )
        unsupervised_training_one_epoch(dataset)

    def test_ZIDatasetCorr_populate_and_train_one(self):
        dataset = ZISyntheticDatasetCorr(n_cells_cluster=10)
        unsupervised_training_one_epoch(dataset)

    def test_corr_zeros(self):
        # Test hierarchy of zeros
        nb_data = SyntheticDatasetCorr()
        zi_data = ZISyntheticDatasetCorr()
        zi_zeros_frac = (zi_data.X == 0).mean()
        nb_zeros_frac = (nb_data.X == 0).mean()

        # nb is not zero inflated
        # zi is zero inflated for all genes
        # We expect the number of zeros to organize accordingly
        self.assertLess(nb_zeros_frac, zi_zeros_frac)
        # We enforce that the zero inflated model has at least 20% of zeros
        self.assertGreaterEqual(zi_zeros_frac, 0.2)
