from unittest import TestCase

from scvi.dataset import SmfishDataset
from . import unsupervised_training_one_epoch


class TestSmfishDataset(TestCase):
    def test_populate(self):
        dataset = SmfishDataset(cell_type_level="minor")
        self.assertEqual(dataset.cell_types[0], "Excluded")
        self.assertEqual(dataset.cell_types[1], "Pyramidal L6")

    def test_train_one(self):
        dataset = SmfishDataset(cell_type_level="minor")
        unsupervised_training_one_epoch(dataset)
