from unittest import TestCase

from scvi.dataset import BrainLargeDataset
from .utils import unsupervised_training_one_epoch


class TestBrainLargeDataset(TestCase):
    def test_populate(self):
        dataset = BrainLargeDataset(
            save_path="tests/data",
            sample_size_gene_var=10,
            nb_genes_to_keep=10,
            max_cells_to_keep=128
        )
        unsupervised_training_one_epoch(dataset)
