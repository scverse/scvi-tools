from unittest import TestCase

import numpy as np

from scvi.dataset import GeneExpressionDataset


class TestGeneExpressionDataset(TestCase):
    def test_X(self):
        data = np.ones((25, 10)) * 100
        dataset = GeneExpressionDataset(data)

        self.assertEqual(dataset.nb_genes, 10)
        self.assertEqual(dataset.nb_cells, 25)

    def test_populate_from_batch_array(self):
        data = np.random.rand(3, 7, 10)
        dataset = GeneExpressionDataset()
        dataset.populate_from_per_batch_array(data)
        self.assertEqual(dataset.nb_cells, 75)
        self.assertEqual(dataset.nb_genes, 10)

    def test_populate_from_per_batch_list(self):
        data = [np.random.rand(7, 10), np.random.rand(5, 10), np.random.rand(3, 10)]
        dataset = GeneExpressionDataset()
        dataset.populate_from_per_batch_list(data)
        self.assertEqual(dataset.nb_cells, 15)
        self.assertEqual(dataset.nb_genes, 10)
        true_batch_indices = np.concatenate([np.zeros(7), np.ones(5), 2 * np.ones(3)])
        self.assertListEqual(dataset.batch_indices, true_batch_indices)

    def test_populate_from_per_label_list(self):
        data = [np.random.rand(7, 10), np.random.rand(5, 10), np.random.rand(3, 10)]
        dataset = GeneExpressionDataset()
        dataset.populate_from_per_batch_list(data)
        self.assertEqual(dataset.nb_cells, 15)
        self.assertEqual(dataset.nb_genes, 10)
        true_labels = np.concatenate([np.zeros(7), np.ones(5), 2 * np.ones(3)])
        self.assertListEqual(dataset.labels, true_labels)

    def test_populate_from_datasets(self):
        data = np.random(5, 10)
        gene_names = np.arange(10).astype(str)
        data1 = np.random(7, 3)
        gene_names1 = np.arange(3).astype(str)
        data2 = np.random(2, 5)
        gene_names1 = np.arange(5).astype(str)
        dataset = GeneExpressionDataset()
        dataset.populate_from_datasets(dataset)

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
