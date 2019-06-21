from unittest import TestCase

from scvi.dataset import HematoDataset
from .utils import unsupervised_training_one_epoch


class TestHematoDataset(TestCase):
    def test_populate(self):
        dataset = HematoDataset(save_path="data/tests/HEMATO")
        unsupervised_training_one_epoch(dataset)
