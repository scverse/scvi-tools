from unittest import TestCase

from scvi.dataset import CortexDataset
from .utils import unsupervised_training_one_epoch


class TestCortexDataset(TestCase):
    def test_populate(self):
        dataset = CortexDataset(save_path="tests/data")
        unsupervised_training_one_epoch(dataset)
