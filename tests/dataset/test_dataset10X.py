from unittest import TestCase

from scvi.dataset.dataset10X import Dataset10X, BrainSmallDataset, available_datasets
from .utils import unsupervised_training_one_epoch


class TestDataset10X(TestCase):
    def test_populate_and_train_one_v1(self):
        dataset = Dataset10X(
            dataset_name=available_datasets["1.1.0"][0], remove_extracted_data=True
        )
        unsupervised_training_one_epoch(dataset)

    def test_brain_small(self):
        dataset = BrainSmallDataset(
            save_path="tests/data",
            save_path_10X="tests/dataset/10X",
            remove_extracted_data=True,
        )
        unsupervised_training_one_epoch(dataset)
