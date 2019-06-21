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

    def test_populate_from_batch_array(self):
        data = np.random.randint(1, 5, size=(3, 7, 10))
        dataset = GeneExpressionDataset()
        dataset.populate_from_per_batch_array(data)
        self.assertEqual(dataset.nb_cells, 75)
        self.assertEqual(dataset.nb_genes, 10)

    def test_populate_from_per_batch_list(self):
        data = [
            np.random.randint(1, 5, size=(7, 10)),
            np.random.randint(1, 5, size=(5, 10)),
            np.random.randint(1, 5, size=(3, 10)),
        ]
        dataset = GeneExpressionDataset()
        dataset.populate_from_per_batch_list(data)
        self.assertEqual(dataset.nb_cells, 15)
        self.assertEqual(dataset.nb_genes, 10)
        true_batch_indices = np.concatenate(
            [
                np.zeros((7, 1), dtype=int),
                np.ones((5, 1), dtype=int),
                2 * np.ones((3, 1), dtype=int),
            ]
        )
        self.assertListEqual(
            dataset.batch_indices.tolist(), true_batch_indices.tolist()
        )

    def test_populate_from_per_label_list(self):
        data = [
            np.random.randint(1, 5, size=(7, 10)),
            np.random.randint(1, 5, size=(5, 10)),
            np.random.randint(1, 5, size=(3, 10)),
        ]
        dataset = GeneExpressionDataset()
        dataset.populate_from_per_label_list(data)
        self.assertEqual(dataset.nb_cells, 15)
        self.assertEqual(dataset.nb_genes, 10)
        true_labels = np.concatenate(
            [
                np.zeros((7, 1), dtype=int),
                np.ones((5, 1), dtype=int),
                2 * np.ones((3, 1), dtype=int),
            ]
        )
        self.assertListEqual(dataset.labels.tolist(), true_labels.tolist())

    def test_populate_from_datasets(self):
        data1 = np.random.randint(1, 5, size=(5, 10))
        gene_names1 = np.arange(10).astype(str)
        dataset1 = GeneExpressionDataset()
        dataset1.populate_from_data(data1, gene_names=gene_names1)
        data2 = np.random.randint(1, 5, size=(7, 3))
        gene_names2 = np.arange(3).astype(str)
        dataset2 = GeneExpressionDataset()
        dataset2.populate_from_data(data2, gene_names=gene_names2)
        data3 = np.random.randint(1, 5, size=(2, 5))
        gene_names3 = np.arange(5).astype(str)
        dataset3 = GeneExpressionDataset()
        dataset3.populate_from_data(data3, gene_names=gene_names3)

        dataset = GeneExpressionDataset()
        dataset.populate_from_datasets([dataset1, dataset2, dataset3])
        self.assertEqual(dataset.nb_cells, 14)
        self.assertEqual(dataset.nb_genes, 3)
        self.assertListEqual(dataset.gene_names.tolist(), ["0", "1", "2"])
        # test for cell_types handling

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
        self.assertEquals(dataset.gene_names[0], "gene_98")

    def test_filter_genes(self):
        pass

    def test_normalize(self):
        pass

    def test_corrupt(self):
        pass
