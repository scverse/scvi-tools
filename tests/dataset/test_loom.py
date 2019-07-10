from unittest import TestCase

from scvi.dataset import LoomDataset, RetinaDataset
from .utils import unsupervised_training_one_epoch


class TestLoomDataset(TestCase):
    def test_populate(self):
        dataset = LoomDataset(filename="Cortex.loom", save_path="tests/data")
        unsupervised_training_one_epoch(dataset)


class TestRetinaDataset(TestCase):
    def test_train_one(self):
        dataset = RetinaDataset(save_path="tests/data")
        unsupervised_training_one_epoch(dataset)
