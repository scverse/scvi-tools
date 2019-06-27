from unittest import TestCase

import anndata
import numpy as np

from scvi.dataset import AnnDatasetFromAnnData, DownloadableAnnDataset
from . import unsupervised_training_one_epoch


class TestAnnDatasetFromAnnData(TestCase):
    def test_init(self):
        data = np.random.randint(1, 5, size=(3, 7))
        ad = anndata.AnnData(data)
        dataset = AnnDatasetFromAnnData(ad)
        self.assertEqual(3, dataset.nb_cells)
        self.assertEqual(7, dataset.nb_genes)

    def test_train_one(self):
        data = np.random.randint(1, 5, size=(4, 7))
        ad = anndata.AnnData(data)
        dataset = AnnDatasetFromAnnData(ad)
        unsupervised_training_one_epoch(dataset)


class TestDownloadableAnnDataset(TestCase):
    def test_populate_and_train_one(self):
        dataset = DownloadableAnnDataset("TM_droplet_mat.h5ad", save_path="tests/data")
        unsupervised_training_one_epoch(dataset)
