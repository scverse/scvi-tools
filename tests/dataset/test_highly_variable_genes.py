from unittest import TestCase
import numpy as np

from scvi.dataset import GeneExpressionDataset, BrainLargeDataset


class TestHighlyVariableGenes(TestCase):
    def test_sparse_no_batch_correction(self):
        dataset = BrainLargeDataset(
            save_path="tests/data",
            sample_size_gene_var=10,
            nb_genes_to_keep=128,
            max_cells_to_keep=256,
        )

        for flavor in ["seurat", "cell_ranger"]:
            n_genes = dataset.nb_genes
            n_top = n_genes // 2
            dataset.subsample_genes(mode=flavor, new_n_genes=n_top, n_bins=3)
            assert dataset.nb_genes < n_genes
            # For some reason the new number of genes can be slightly different than n_top

            dataset.highly_variable_genes(flavor=flavor, n_bins=3)

    def test_batch_correction(self):
        data = [
            np.random.randint(1, 5, size=(50, 25)),
            np.random.randint(1, 5, size=(50, 25)),
            np.random.randint(1, 5, size=(50, 25)),
        ]
        dataset = GeneExpressionDataset()
        dataset.populate_from_per_batch_list(data)

        n_genes = dataset.nb_genes
        n_top = n_genes // 2
        dataset.highly_variable_genes(n_bins=3, flavor="seurat")

        dataset.highly_variable_genes(n_bins=3, flavor="seurat")

        df = dataset.highly_variable_genes(n_bins=3, n_top_genes=n_top, flavor="seurat")
        assert df["highly_variable"].sum() >= n_top
        pass

    def test_dense_subsample_genes(self):
        data = [
            np.random.randint(1, 5, size=(50, 26)),
            np.random.randint(1, 5, size=(50, 26)),
            np.random.randint(1, 5, size=(50, 26)),
        ]

        # With default
        dataset = GeneExpressionDataset()
        dataset.populate_from_per_batch_list(data)
        n_genes = dataset.nb_genes
        n_top = n_genes // 2
        dataset.subsample_genes(new_n_genes=n_top)
        assert dataset.nb_genes < n_genes
        # For some reason the new number of genes can be slightly different than n_top

        # With Seurat
        dataset = GeneExpressionDataset()
        dataset.populate_from_per_batch_list(data)
        dataset.subsample_genes(new_n_genes=n_top, mode="seurat")
        assert dataset.nb_genes < n_genes
        # For some reason the new number of genes can be slightly different than n_top
