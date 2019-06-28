from unittest import TestCase

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


class TestSyntheticRandomDataset(TestCase):
    def test_populate_and_train_one(self):
        dataset = SyntheticRandomDataset(save_path="tests/data")
        unsupervised_training_one_epoch(dataset)


class TestSyntheticDatasetCorr(TestCase):
    def test_populate_and_train_one(self):
        dataset = SyntheticDatasetCorr(n_cells_cluster=10)
        unsupervised_training_one_epoch(dataset)


class TestZISyntheticDatasetCorr(TestCase):
    def test_populate_and_train_one(self):
        dataset = ZISyntheticDatasetCorr(n_cells_cluster=10)
        unsupervised_training_one_epoch(dataset)
