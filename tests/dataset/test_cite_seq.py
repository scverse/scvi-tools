from unittest import TestCase

from scvi.dataset import CiteSeqDataset
from .utils import unsupervised_training_one_epoch


class TestCiteSeqDataset(TestCase):
    def test_populate_and_train_one(self):
        for name in ["cbmc", "pbmc"]:
            dataset = CiteSeqDataset(name=name, save_path="tests/data/citeSeq/")
            unsupervised_training_one_epoch(dataset)
