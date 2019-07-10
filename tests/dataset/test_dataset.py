from unittest import TestCase

import numpy as np

from scvi.dataset import CortexDataset, GeneExpressionDataset
from scvi.dataset.dataset import remap_categories, CellMeasurement


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

    def test_populate_from_data_with_measurements(self):
        data = np.ones((25, 10)) * 100
        paired = np.ones((25, 4)) * np.arange(0, 4)
        pair_names = ["gabou", "achille", "pedro", "oclivio"]
        y = CellMeasurement(
            name="dev", data=paired, columns_attr_name="dev_names", columns=pair_names
        )
        dataset = GeneExpressionDataset()

        dataset.populate_from_data(data, Ys=[y])

        self.assertEqual(dataset.nb_genes, 10)
        self.assertEqual(dataset.nb_cells, 25)

        self.assertTrue(hasattr(dataset, "dev"))
        self.assertTrue(hasattr(dataset, "dev_names"))

        self.assertListEqual(dataset.dev_names.tolist(), pair_names)
        self.assertListEqual(dataset.dev[0].tolist(), [0, 1, 2, 3])

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
            ["GENE_0", "GENE_1", "GENE_2"], dataset.gene_names.tolist()
        )

        # test for labels sharing
        dataset2.labels = [0, 0, 0, 1, 1, 1, 1]
        dataset2.initialize_mapped_attribute("labels", "cell_types", ["0", "1"])
        dataset3.labels = [0, 1]
        dataset3.initialize_mapped_attribute("labels", "cell_types", ["0", "2"])
        dataset = GeneExpressionDataset()
        dataset.populate_from_datasets([dataset2, dataset3], shared_labels=True)
        self.assertListEqual(
            np.squeeze(dataset.labels).tolist(), [0, 0, 0, 1, 1, 1, 1, 0, 2]
        )
        self.assertListEqual(dataset.cell_types, ["0", "1", "2"])

        dataset_unshared = GeneExpressionDataset()
        dataset_unshared.populate_from_datasets(
            [dataset2, dataset3], shared_labels=False
        )
        self.assertListEqual(
            np.squeeze(dataset_unshared.labels).tolist(), [0, 0, 0, 1, 1, 1, 1, 2, 3]
        )
        self.assertListEqual(dataset_unshared.cell_types, ["0", "1", "0", "2"])

        # test for batch_indices offsetting
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
            np.squeeze(dataset.batch_indices).tolist(), [0, 0, 0, 1, 1, 1, 1, 2, 3]
        )
        self.assertListEqual(
            getattr(dataset, "experiment"), ["fish_2", "scrna_2", "fish_3", "scrna_3"]
        )

    def test_populate_from_datasets_gene_attributes_merging(self):
        data = np.random.randint(1, 5, size=(5, 10))
        gene_names = np.array(["gene_%d" % i for i in range(10)])
        gene_attr1 = np.array([["1"] for _ in range(10)])
        gene_attr2 = np.array([["2"] for _ in range(10)])
        dataset1 = GeneExpressionDataset()
        dataset2 = GeneExpressionDataset()

        dataset1.populate_from_data(
            data, gene_names=gene_names, gene_attributes_dict={"test": gene_attr1}
        )
        dataset2.populate_from_data(
            data, gene_names=gene_names, gene_attributes_dict={"test": gene_attr2}
        )

        dataset = GeneExpressionDataset()
        dataset.populate_from_datasets([dataset1, dataset2])

        # Should keep the gene attribute of the first dataset
        self.assertEqual(dataset.test[0, 0], "1")

    def test_populate_from_datasets_cell_attributes_merging(self):
        data = np.random.randint(1, 5, size=(5, 10))
        gene_names = np.array(["gene_%d" % i for i in range(10)])
        cell_attr1 = np.array([["1"] for _ in range(5)])
        cell_attr2 = np.array([["2"] for _ in range(5)])
        dataset1 = GeneExpressionDataset()
        dataset2 = GeneExpressionDataset()

        dataset1.populate_from_data(
            data, gene_names=gene_names, cell_attributes_dict={"test": cell_attr1}
        )
        dataset2.populate_from_data(
            data, gene_names=gene_names, cell_attributes_dict={"test": cell_attr2}
        )

        dataset = GeneExpressionDataset()
        dataset.populate_from_datasets([dataset1, dataset2])
        self.assertTupleEqual(dataset.test.shape, (10, 1))
        self.assertListEqual(np.squeeze(dataset.test).tolist(), ["1"] * 5 + ["2"] * 5)

    def test_populate_from_datasets_with_measurments(self):
        data = np.random.randint(1, 5, size=(5, 10))
        gene_names = np.array(["gene_%d" % i for i in range(10)])

        paired1 = np.ones((5, 5)) * np.arange(0, 5)
        pair_names1 = ["gabou", "achille", "pedro", "oclivio", "gayoso"]
        y1 = CellMeasurement(
            name="dev", data=paired1, columns_attr_name="dev_names", columns=pair_names1
        )
        paired2 = np.ones((5, 4)) * np.arange(0, 4)
        pair_names2 = ["gabou", "oclivio", "achille", "pedro"]
        y2 = CellMeasurement(
            name="dev", data=paired2, columns_attr_name="dev_names", columns=pair_names2
        )

        dataset1 = GeneExpressionDataset()
        dataset2 = GeneExpressionDataset()

        dataset1.populate_from_data(data, Ys=[y1], gene_names=gene_names)
        dataset2.populate_from_data(data, Ys=[y2], gene_names=gene_names)

        dataset = GeneExpressionDataset()
        dataset.populate_from_datasets([dataset1, dataset2])

        self.assertTrue(hasattr(dataset, "dev"))
        self.assertTrue(hasattr(dataset, "dev_names"))

        self.assertListEqual(
            dataset.dev_names.tolist(), ["achille", "gabou", "oclivio", "pedro"]
        )
        self.assertListEqual(dataset.dev[0].tolist(), [1, 0, 3, 2])
        self.assertListEqual(dataset.dev[5].tolist(), [2, 0, 1, 3])

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
        self.assertTupleEqual((25, 1), dataset.labels.shape)
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
        self.assertEqual(dataset.gene_names[0], "GENE_99")
        dataset.subsample_genes(subset_genes=[1, 6, 7])
        self.assertEqual(dataset.gene_names[0], "GENE_98")

    def test_filter_genes(self):
        data = np.random.randint(1, 5, size=(5, 10))
        gene_names = np.array(["gene_%d" % i for i in range(10)])

        dataset = GeneExpressionDataset()
        dataset.populate_from_data(data, gene_names=gene_names)
        gene_names_true = ["GENE_1", "GENE_3"]
        dataset.filter_genes_by_attribute(gene_names_true)
        self.assertListEqual(gene_names_true, dataset.gene_names.tolist())

    def test_reorder_genes(self):
        data = np.ones((25, 100)) * 100

        gene_names = np.array(["gene_%d" % i for i in range(100)])
        dataset = GeneExpressionDataset()
        dataset.populate_from_data(data, gene_names=gene_names)
        dataset.reorder_genes(["GENE_2", "GENE_47"])
        # New order should be 2, 47, 0, 1, 3
        self.assertListEqual(
            list(dataset.gene_names[0:5]),
            ["GENE_2", "GENE_47", "GENE_0", "GENE_1", "GENE_3"],
        )

    def test_genes_to_idx(self):
        data = np.random.randint(1, 5, size=(5, 10))
        gene_names = np.array(["gene_%d" % i for i in range(10)])

        dataset = GeneExpressionDataset()
        dataset.populate_from_data(data, gene_names=gene_names)
        indices = dataset.genes_to_index(["GENE_%d" % i for i in range(10)])
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
        dataset.filter_cells_by_count()
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
        data = np.random.randint(1, 5, size=(8, 20))
        labels = [0, 0, 1, 2, 2, 1, 0, 1]
        cell_types = ["0", "1", "2"]

        dataset = GeneExpressionDataset()
        dataset.populate_from_data(data, labels=labels, cell_types=cell_types)
        dataset.merge_cell_types(["0", "1"], new_cell_type_name="0 and 1")
        self.assertListEqual(
            [[3], [3], [3], [2], [2], [3], [3], [3]], dataset.labels.tolist()
        )
        dataset.remap_categorical_attributes()
        self.assertListEqual(
            [[1], [1], [1], [0], [0], [1], [1], [1]], dataset.labels.tolist()
        )
        self.assertListEqual(["2", "0 and 1"], dataset.cell_types.tolist())

    def test_map_cell_types(self):
        data = np.random.randint(1, 5, size=(7, 10))
        labels = [0, 0, 4, 4, 2, 3, 5]
        cell_types = ["0", "1", "2", "3", "4", "5"]

        dataset = GeneExpressionDataset()
        dataset.populate_from_data(data, labels=labels, cell_types=cell_types)
        dataset.map_cell_types({("0", "2"): "6", ("3", "4"): "7"})
        dataset.remap_categorical_attributes()
        self.assertListEqual(dataset.cell_types.tolist(), ["5", "6", "7"])
        self.assertListEqual(np.squeeze(dataset.labels).tolist(), [1, 1, 2, 2, 1, 2, 0])


class TestCollate(TestCase):
    def test_collate_normal(self):
        data = np.ones((25, 2)) * np.arange(0, 25).reshape((-1, 1))
        batch_indices = np.arange(0, 25).reshape((-1, 1))
        dataset = GeneExpressionDataset()
        dataset.populate_from_data(data, batch_indices=batch_indices)

        collate_fn = dataset.collate_fn_builder()
        x, mean, var, batch, labels = collate_fn([1, 2])
        self.assertListEqual(x.tolist(), [[1.0, 1.0], [2.0, 2.0]])
        self.assertListEqual(batch.tolist(), [[1], [2]])

    def test_collate_add(self):
        data = np.ones((25, 2)) * np.arange(0, 25).reshape((-1, 1))
        batch_indices = np.arange(0, 25).reshape((-1, 1))
        x_coords = np.arange(0, 25).reshape((-1, 1))
        proteins = (
            np.ones((25, 3)) + np.arange(0, 25).reshape((-1, 1)) + np.arange(0, 3)
        )
        proteins_name = ["A", "B", "C"]
        dataset = GeneExpressionDataset()
        dataset.populate_from_data(
            data,
            batch_indices=batch_indices,
            cell_attributes_dict={"x_coords": x_coords},
            Ys=[
                CellMeasurement(
                    name="proteins",
                    data=proteins,
                    columns_attr_name="protein_names",
                    columns=proteins_name,
                )
            ],
        )

        collate_fn = dataset.collate_fn_builder(
            add_attributes_and_types={"x_coords": np.float32, "proteins": np.float32}
        )
        x, mean, var, batch, labels, x_coords_tensor, proteins_tensor = collate_fn(
            [1, 2]
        )
        self.assertListEqual(x_coords_tensor.tolist(), [[1.0], [2.0]])
        self.assertListEqual(
            proteins_tensor.tolist(), [[2.0, 3.0, 4.0], [3.0, 4.0, 5.0]]
        )


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
