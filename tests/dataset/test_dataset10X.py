from unittest import TestCase

from scvi.dataset.dataset10X import Dataset10X, BrainSmallDataset
from .utils import unsupervised_training_one_epoch


class TestDataset10X(TestCase):
    def test_populate_and_train_one_v1(self):
        dataset = Dataset10X(
            dataset_name="cd4_t_helper",
            remove_extracted_data=True,
            save_path="tests/data/10X",
        )
        unsupervised_training_one_epoch(dataset)

    def test_brain_small(self):
        dataset = BrainSmallDataset(
            save_path="tests/data",
            save_path_10X="tests/data/10X",
            remove_extracted_data=True,
        )
        unsupervised_training_one_epoch(dataset)
