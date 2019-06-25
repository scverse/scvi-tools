from unittest import TestCase

import numpy as np

from scvi.dataset import CortexDataset, GeneExpressionDataset
from scvi.dataset.dataset import remap_categories


class TestGeneExpressionDataset(TestCase):
    def test_populate_from_data(self):
        data = np.ones((25, 10)) * 100
        dataset = GeneExpressionDataset()
        dataset.populate_from_data(data)

        self.assertEqual(dataset.nb_genes, 10)
        self.assertEqual(dataset.nb_cells, 25)
        # default batch_indices and labels
        self.assertListEqual([[0] for i in range(25)], dataset.batch_indices.tolist())
        self.assertListEqual([[0] for i in range(25)], dataset.labels.tolist())

    def test_populate_from_batch_array(self):
        data = np.random.randint(1, 5, size=(3, 7, 10))
        dataset = GeneExpressionDataset()
        dataset.populate_from_per_batch_array(data)
        self.assertEqual(dataset.nb_cells, 21)
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
            true_batch_indices.tolist(), dataset.batch_indices.tolist()
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

    def test_populate_from_datasets_dummy_data(self):
        data1 = np.random.randint(1, 5, size=(5, 10))
        gene_names1 = np.array(["gene_%d" % i for i in range(10)])
        dataset1 = GeneExpressionDataset()
        dataset1.populate_from_data(data1, gene_names=gene_names1)
        data2 = np.random.randint(1, 5, size=(7, 3))
        gene_names2 = np.array(["gene_%d" % i for i in range(3)])
        dataset2 = GeneExpressionDataset()
        dataset2.populate_from_data(data2, gene_names=gene_names2)
        data3 = np.random.randint(1, 5, size=(2, 5))
        gene_names3 = np.array(["gene_%d" % i for i in range(5)])
        dataset3 = GeneExpressionDataset()
        dataset3.populate_from_data(data3, gene_names=gene_names3)

        dataset = GeneExpressionDataset()
        dataset.populate_from_datasets([dataset1, dataset2, dataset3])
        self.assertEqual(14, dataset.nb_cells)
        self.assertEqual(3, dataset.nb_genes)
        self.assertListEqual(
            ["gene_0", "gene_1", "gene_2"], dataset.gene_names.tolist()
        )

        # test for attribute mapping handling
        # sharing instruction "concatenation"
        dataset2.labels = [0, 0, 0, 1, 1, 1, 1]
        dataset2.initialize_mapped_attribute("labels", "cell_types", ["0", "1"])
        dataset3.labels = [0, 1]
        dataset3.initialize_mapped_attribute("labels", "cell_types", ["0", "2"])
        dataset = GeneExpressionDataset()
        dataset.populate_from_datasets([dataset2, dataset3])
        self.assertListEqual(
            [0, 0, 0, 1, 1, 1, 1, 0, 2], np.squeeze(dataset.labels).tolist()
        )
        self.assertListEqual(["0", "1", "2"], dataset.cell_types)

        # sharing instruction "offset"
        dataset2.batch_indices = [0, 0, 0, 1, 1, 1, 1]
        dataset2.initialize_mapped_attribute(
            "batch_indices", "experiment", ["fish_2", "scrna_2"]
        )
        dataset3.batch_indices = [0, 1]
        dataset3.initialize_mapped_attribute(
            "batch_indices", "experiment", ["fish_3", "scrna_3"]
        )
        dataset = GeneExpressionDataset()
        dataset.populate_from_datasets([dataset2, dataset3])
        self.assertListEqual(
            [0, 0, 0, 1, 1, 1, 1, 2, 3], np.squeeze(dataset.batch_indices).tolist()
        )
        self.assertListEqual(
            ["fish_2", "scrna_2", "fish_3", "scrna_3"], getattr(dataset, "experiment")
        )

    def test_populate_from_datasets_cortex(self):
        cortex_dataset_1 = CortexDataset(save_path="tests/data")
        cortex_dataset_1.subsample_genes(subset_genes=np.arange(0, 3))
        cortex_dataset_1.filter_cell_types(["microglia", "oligodendrocytes"])
        cortex_dataset_2 = CortexDataset(save_path="tests/data")
        cortex_dataset_2.subsample_genes(subset_genes=np.arange(1, 4))
        cortex_dataset_2.filter_cell_types(
            ["endothelial-mural", "interneurons", "microglia", "oligodendrocytes"]
        )
        cortex_dataset_2.filter_cell_types([2, 0])
        dataset = GeneExpressionDataset()
        dataset.populate_from_datasets([cortex_dataset_1, cortex_dataset_2])
        self.assertEqual(2, dataset.nb_genes)

    def test_labels(self):
        data = np.ones((25, 10)) * 100
        labels = np.array(range(25))
        dataset = GeneExpressionDataset()
        dataset.populate_from_data(data, labels=labels)
        self.assertTupleEqual((25, 1), dataset.labels.shape, )
        self.assertEqual(dataset.labels[5, 0], 5)

        labels = np.ones(25) * 5
        dataset = GeneExpressionDataset()
        dataset.populate_from_data(data, labels=labels)
        self.assertTupleEqual(dataset.labels.shape, (25, 1))
        self.assertEqual(dataset.labels[5, 0], 0)

    def test_remap_categorical_attributes(self):
        data = np.random.randint(1, 5, size=(7, 11))
        labels = [1, 1, 1, 1, 1, 2, 2]
        dataset = GeneExpressionDataset()
        dataset.populate_from_data(data, labels=labels)

        labels_true = [0, 0, 0, 0, 0, 1, 1]
        labels_true = [[i] for i in labels_true]
        self.assertListEqual(labels_true, dataset.labels.tolist())

    def test_compute_library_size_batch(self):
        data = np.exp(10) / 10 * np.ones((7, 10), dtype=int)
        dataset = GeneExpressionDataset()
        dataset.populate_from_data(data)

        local_means_true = [[10.0] for _ in range(7)]
        local_vars_true = [[0.0] for _ in range(7)]
        self.assertEqual(local_means_true, dataset.local_means.tolist())
        self.assertEqual(local_vars_true, dataset.local_vars.tolist())

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
        # The most variable genes should be in first position
        self.assertEquals(dataset.gene_names[0], "gene_99")
        dataset.subsample_genes(subset_genes=[1, 6, 7])
        self.assertEquals(dataset.gene_names[0], "gene_98")

    def test_filter_genes(self):
        data = np.random.randint(1, 5, size=(5, 10))
        gene_names = np.array(["gene_%d" % i for i in range(10)])

        dataset = GeneExpressionDataset()
        dataset.populate_from_data(data, gene_names=gene_names)
        gene_names_true = ["gene_1", "gene_3"]
        dataset.filter_genes(gene_names_true)
        self.assertListEqual(gene_names_true, dataset.gene_names.tolist())

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

    def test_genes_to_idx(self):
        data = np.random.randint(1, 5, size=(5, 10))
        gene_names = np.array(["gene_%d" % i for i in range(10)])

        dataset = GeneExpressionDataset()
        dataset.populate_from_data(data, gene_names=gene_names)
        indices = dataset.genes_as_index(["gene_%d" % i for i in range(10)])
        self.assertListEqual([i for i in range(10)], indices.tolist())

    def test_subsample_cells(self):
        data = np.arange(1, 6)[:, None] * np.ones(7)[None, :]

        dataset = GeneExpressionDataset()
        dataset.populate_from_data(data)
        # default
        dataset.subsample_cells()
        self.assertEqual(5, dataset.nb_cells)
        # when size is a float
        dataset.subsample_cells(size=0.8)
        data_true = np.arange(5, 1, -1)[:, None] * np.ones(7)[None, :]
        self.assertListEqual(data_true.tolist(), dataset.X.tolist())
        # when size is an int
        dataset.subsample_cells(size=2)
        self.assertEqual(2, dataset.nb_cells)

    def test_filter_cells(self):
        data = np.ones((25, 10)) * 100
        data[4:6, :] = 0
        dataset = GeneExpressionDataset()
        dataset.populate_from_data(data)
        self.assertEqual(25, dataset.nb_cells)
        dataset.filter_cells()
        self.assertEqual(23, dataset.nb_cells)

    def test_filter_cell_types(self):
        data = np.random.randint(1, 5, size=(5, 10))
        labels = [0, 0, 1, 1, 1]
        cell_types = ["0", "1"]

        dataset = GeneExpressionDataset()
        dataset.populate_from_data(data, labels=labels, cell_types=cell_types)
        dataset.filter_cell_types(["0"])
        self.assertListEqual(data[:2].tolist(), dataset.X.tolist())

    def test_merge_cell_types(self):
        data = np.random.randint(1, 5, size=(5, 10))
        labels = [0, 0, 1, 1, 1]
        cell_types = ["0", "1"]

        dataset = GeneExpressionDataset()
        dataset.populate_from_data(data, labels=labels, cell_types=cell_types)
        dataset.merge_cell_types(["0", "1"], new_cell_type_name="0 and 1")
        dataset.remap_categorical_attributes()
        self.assertListEqual([[0] for _ in range(5)], dataset.labels.tolist())
        self.assertListEqual(["0 and 1"], dataset.cell_types.tolist())

    def test_map_cell_types(self):
        data = np.random.randint(1, 5, size=(5, 10))
        labels = [0, 0, 1, 1, 2]
        cell_types = ["0", "1", "2"]

        dataset = GeneExpressionDataset()
        dataset.populate_from_data(data, labels=labels, cell_types=cell_types)
        dataset.map_cell_types({("0", "2"): "3"})


class TestRemapCategories(TestCase):
    def test_remap_categories(self):
        labels = [0, 0, 0, 2, 2, 3]
        labels, n_labels = remap_categories(labels)
        labels_true = [0, 0, 0, 1, 1, 2]
        self.assertListEqual(labels_true, labels.tolist())
        self.assertEqual(3, n_labels)

        # with absent categories and mappings
        labels = [2, 2, 3]
        mappings_dict = {"cell_types": ["0", "1", "2", "3"]}
        labels, n_labels, mappings = remap_categories(
            labels, mappings_dict=mappings_dict
        )
        labels_true = [0, 0, 1]
        self.assertListEqual(labels_true, labels.tolist())
        self.assertEqual(2, n_labels)
        self.assertListEqual(["2", "3"], mappings["cell_types"].tolist())


class TestSpatialGeneExpressionDataset(TestCase):
    def test_make_tensor_batch_from_indices(self):
        pass
