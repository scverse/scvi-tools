from unittest import TestCase

from scvi.dataset import SeqfishDataset
from .utils import unsupervised_training_one_epoch


class TestSeqfishDataset(TestCase):
    def test_populate(self):
        dataset = SeqfishDataset(save_path="tests/data")
        unsupervised_training_one_epoch(dataset)
