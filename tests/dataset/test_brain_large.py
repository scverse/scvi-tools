from unittest import TestCase

from scvi.dataset import BrainLargeDataset
from . import unsupervised_training_one_epoch


class TestBrainLargeDataset(TestCase):
    def test_populate(self):
        dataset = BrainLargeDataset(save_path="tests/data")
        unsupervised_training_one_epoch(dataset)
