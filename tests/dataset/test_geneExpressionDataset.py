from unittest import TestCase

import numpy as np

from scvi.dataset import GeneExpressionDataset


class TestGeneExpressionDataset(TestCase):
    def test_X(self):
        data = np.ones((25, 10)) * 100
        dataset = GeneExpressionDataset()
        dataset.populate_from_data(data)

        self.assertEqual(dataset.nb_genes, 10)
        self.assertEqual(dataset.nb_cells, 25)

    def test_filter_cells(self):
        data = np.ones((25, 10)) * 100
        data[4:6, :] = 0
        dataset = GeneExpressionDataset()
        dataset.populate_from_data(data)
        self.assertEqual(dataset.nb_cells, 25)
        dataset.filter_cells()
        self.assertEqual(dataset.nb_cells, 23)

    def test_labels(self):
        data = np.ones((25, 10)) * 100
        labels = np.array(range(25))
        dataset = GeneExpressionDataset()
        dataset.populate_from_data(data, labels=labels)
        self.assertTupleEqual(dataset.labels.shape, (25, 1))
        self.assertEqual(dataset.labels[5, 0], 5)

        labels = np.ones(25) * 5
        dataset = GeneExpressionDataset()
        dataset.populate_from_data(data, labels=labels)
        self.assertTupleEqual(dataset.labels.shape, (25, 1))
        self.assertEqual(dataset.labels[5, 0], 0)

    def test_update_genes(self):
        pass

    def test_update_cells(self):
        pass

    def test_subsample_genes(self):
        data = np.ones((25, 100)) * 100
        variable_data = data
        variable_data[0, :] = 2
        variable_data *= np.arange(0, 100)

        gene_names = np.array(["gene_%d" % i for i in range(100)])
        dataset = GeneExpressionDataset()
        dataset.populate_from_data(data, gene_names=gene_names)
        dataset.subsample_genes(new_ratio_genes=0.4)
        self.assertTupleEqual(dataset.gene_names.shape, (40,))
        dataset.subsample_genes(new_n_genes=25)
        self.assertTupleEqual(dataset.gene_names.shape, (25,))
        # Most variable gene should be in first position
        self.assertEquals(dataset.gene_names[0], "gene_99")
        dataset.subsample_genes(subset_genes=[1, 6, 7])
        self.assertEquals(dataset.gene_names[0], "gene_98")

    def test_subsample_cells(self):
        pass

    def test_reorder_genes(self):
        data = np.ones((25, 100)) * 100

        gene_names = np.array(["gene_%d" % i for i in range(100)])
        dataset = GeneExpressionDataset()
        dataset.populate_from_data(data, gene_names=gene_names)
        dataset.reorder_genes(["gene_2", "gene_47"])
        # New order should be 2, 47, 0, 1, 3
        self.assertListEqual(
            list(dataset.gene_names[0:5]),
            ["gene_2", "gene_47", "gene_0", "gene_1", "gene_3"],
        )

    def test_filter_genes(self):
        pass

    def test_normalize(self):
        pass

    def test_corrupt(self):
        pass

    def test_from_batch_array(self):
        pass

    def test_from_per_batch_list(self):
        pass

    def test_from_per_label_list(self):
        pass

    def test_from_datasets(self):
        pass
