from unittest import TestCase

from scvi.dataset import BiomDataset
from .utils import unsupervised_training_one_epoch


class TestBiomDataset(TestCase):
    # fix once NB is fixed
    def test_populate(self):
        dataset = BiomDataset(filename="feature-table.biom",
                              save_path="tests/data")
        unsupervised_training_one_epoch(dataset)

