from unittest import TestCase

from scvi.dataset import CsvDataset
from . import unsupervised_training_one_epoch


class TestCsvDataset(TestCase):
    def test_populate_and_train(self):
        dataset = CsvDataset(
            "GSE100866_CBMC_8K_13AB_10X-RNA_umi.csv.gz",
            save_path="tests/data",
            compression="gzip",
        )
        unsupervised_training_one_epoch(csv_dataset)
