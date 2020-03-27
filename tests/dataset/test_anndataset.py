from unittest import TestCase

import anndata
import numpy as np

from scvi.dataset import (
    AnnDatasetFromAnnData,
    DownloadableAnnDataset,
    CellMeasurement,
    GeneExpressionDataset,
)
from .utils import unsupervised_training_one_epoch
import scipy.sparse as sp_sparse


class TestAnnDataset(TestCase):
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

    def test_populate_and_train_one(self):
        dataset = DownloadableAnnDataset("TM_droplet_mat.h5ad", save_path="tests/data")
        unsupervised_training_one_epoch(dataset)

    def test_use_raw_flag(self):
        raw_data = np.random.randint(1, 5, size=(4, 7))
        ad = anndata.AnnData(raw_data)
        ad.raw = ad.copy()
        dataset = AnnDatasetFromAnnData(ad, use_raw=True)
        np.testing.assert_array_equal(dataset.X, raw_data)

    def test_data_loader(self):
        data = np.ones((25, 10)) * 100
        paired = np.ones((25, 4)) * np.arange(0, 4)
        pair_names = ["gabou", "achille", "pedro", "oclivio"]
        y = CellMeasurement(
            name="dev", data=paired, columns_attr_name="dev_names", columns=pair_names
        )
        dataset = GeneExpressionDataset()
        dataset.populate_from_data(data, Ys=[y])
        ad = dataset.to_anndata()
        dataset_ad = AnnDatasetFromAnnData(
            ad, cell_measurements_col_mappings={"dev": "dev_names"}
        )
        self.assertTrue((paired == dataset_ad.dev).all())
        self.assertTrue((dataset.X == dataset_ad.X).all())
        self.assertTrue((dataset.cell_types == dataset_ad.cell_types).all())

    def test_sparse_data(self):
        data = np.random.poisson(0.2, size=(25, 10))

        sparse_mat = sp_sparse.csr_matrix(data)
        ad = anndata.AnnData(sparse_mat)
        AnnDatasetFromAnnData(ad)

        sparse_mat = sp_sparse.csc_matrix(data)
        ad = anndata.AnnData(sparse_mat)
        AnnDatasetFromAnnData(ad)
