from unittest import TestCase

from scvi.dataset import SmfishDataset
from .utils import unsupervised_training_one_epoch


class TestSmfishDataset(TestCase):
    def test_populate(self):
        dataset = SmfishDataset(use_high_level_cluster=False)
        self.assertEqual(dataset.cell_types[0], "Excluded")
        self.assertEqual(dataset.cell_types[1], "Pyramidal L6")

    def test_train_one(self):
        dataset = SmfishDataset(use_high_level_cluster=False)
        unsupervised_training_one_epoch(dataset)
