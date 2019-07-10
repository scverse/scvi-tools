from unittest import TestCase

from scvi.dataset import CsvDataset, BreastCancerDataset, MouseOBDataset
from .utils import unsupervised_training_one_epoch


class TestCsvDataset(TestCase):
    def test_populate_and_train(self):
        datasets = [
            CsvDataset(
                "GSE100866_CBMC_8K_13AB_10X-RNA_umi.csv.gz",
                save_path="tests/data",
                compression="gzip",
            ),
            BreastCancerDataset(save_path="tests/data"),
            MouseOBDataset(save_path="tests/data"),
        ]
        for dataset in datasets:
            unsupervised_training_one_epoch(dataset)
