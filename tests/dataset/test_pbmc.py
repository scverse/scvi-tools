from unittest import TestCase

from scvi.dataset import PbmcDataset, PurifiedPBMCDataset
from .utils import unsupervised_training_one_epoch


class TestPbmcDataset(TestCase):
    def test_populate(self):
        dataset = PbmcDataset(
            save_path="tests/data/",
            save_path_10X="tests/data/10X",
            remove_extracted_data=True,
        )
        unsupervised_training_one_epoch(dataset)


class TestPurifiedPBMCDataset(TestCase):
    def test_populate(self):
        dataset = PurifiedPBMCDataset(
            save_path="tests/data/10X", remove_extracted_data=True
        )
        unsupervised_training_one_epoch(dataset)
        self.assertEqual(len(dataset.cell_types), 10)

        dataset = PurifiedPBMCDataset(
            save_path="tests/data/10X",
            subset_datasets=list(range(6)),
            remove_extracted_data=True,
        )
        unsupervised_training_one_epoch(dataset)
        self.assertEqual(len(dataset.cell_types), 6)
