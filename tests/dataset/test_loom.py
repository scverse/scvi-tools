from unittest import TestCase

from scvi.dataset import (
    LoomDataset,
    RetinaDataset,
    PreFrontalCortexStarmapDataset,
    FrontalCortexDropseqDataset,
)
from .utils import unsupervised_training_one_epoch


class TestLoomDataset(TestCase):
    def test_populate(self):
        dataset = LoomDataset(filename="Cortex.loom", save_path="tests/data")
        unsupervised_training_one_epoch(dataset)


class TestSubDataset(TestCase):
    def test_retina_load_train_one(self):
        dataset = RetinaDataset(save_path="tests/data")
        unsupervised_training_one_epoch(dataset)

    def test_pfc_starmap_load_train_one(self):
        gene_dataset = PreFrontalCortexStarmapDataset(save_path="tests/data")
        unsupervised_training_one_epoch(gene_dataset)

    def test_fc_dropseq_load_train_one(self):
        gene_dataset = FrontalCortexDropseqDataset(save_path="tests/data")
        unsupervised_training_one_epoch(gene_dataset)
