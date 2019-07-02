from unittest import TestCase

from scvi.dataset import PbmcDataset, PurifiedPBMCDataset
from .utils import unsupervised_training_one_epoch


class TestPbmcDataset(TestCase):
    def test_populate(self):
        dataset = PbmcDataset(save_path="tests/data/", save_path_10X="tests/data/10X")
        unsupervised_training_one_epoch(dataset)


class TestPurifiedPBMCDataset(TestCase):
    def test_populate(self):
        dataset = PurifiedPBMCDataset(save_path="tests/data/10X")
        unsupervised_training_one_epoch(dataset)
