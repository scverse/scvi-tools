from unittest import TestCase

import numpy as np

from scvi.dataset import GeneExpressionDataset


class TestGeneExpressionDataset(TestCase):
    def test_X(self):
        data = np.ones((25, 10)) * 100
        dataset = GeneExpressionDataset(data)

        self.assertEqual(dataset.nb_genes, 10)
        self.assertEqual(dataset.nb_cells, 25)

    def test_filter_cells(self):
        data = np.ones((25, 10)) * 100
        data[4:6, :] = 0
        dataset = GeneExpressionDataset(data)
        self.assertEqual(dataset.nb_cells, 25)
        dataset.filter_cells()
        self.assertEqual(dataset.nb_cells, 23)

    def test_labels(self):
        data = np.ones((25, 10)) * 100
        labels = np.array(range(25))
        dataset = GeneExpressionDataset(data, labels=labels)
        self.assertTupleEqual(dataset.labels.shape, (25, 1))
        self.assertEqual(dataset.labels[5, 0], 5)

        labels = np.ones(25) * 5
        dataset = GeneExpressionDataset(data, labels=labels)
        self.assertTupleEqual(dataset.labels.shape, (25, 1))
        self.assertEqual(dataset.labels[5, 0], 0)

    def test_update_genes(self):
        self.fail()

    def test_update_cells(self):
        self.fail()

    def test_subsample_genes(self):
        self.fail()

    def test_subsample_cells(self):
        self.fail()

    def test_reorder_genes(self):
        self.fail()

    def test_filter_genes(self):
        self.fail()

    def test_normalize(self):
        self.fail()

    def test_corrupt(self):
        self.fail()

    def test_from_batch_array(self):
        self.fail()

    def test_from_per_batch_list(self):
        self.fail()

    def test_from_per_label_list(self):
        self.fail()

    def test_from_datasets(self):
        self.fail()
