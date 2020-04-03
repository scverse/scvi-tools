from unittest import TestCase
import numpy as np
from scvi.dataset import BrainLargeDataset, SyntheticDataset


class TestHighlyVariableGenes(TestCase):
    def test_sparse_no_batch_correction(self):
        for flavor in ["seurat_v2", "cell_ranger", "seurat_v3"]:
            dataset = BrainLargeDataset(
                save_path="tests/data",
                sample_size_gene_var=10,
                nb_genes_to_keep=128,
                max_cells_to_keep=256,
            )

            n_genes = dataset.nb_genes
            n_top = n_genes // 2
            dataset.subsample_genes(mode=flavor, new_n_genes=n_top, n_bins=3)
            assert dataset.nb_genes < n_genes
            # For some reason the new number of genes can be slightly different than n_top

            dataset._highly_variable_genes(
                flavor=flavor,
                n_bins=3,
                n_top_genes=3 if flavor == "seurat_v3" else None,
            )

    def test_batch_correction(self):
        dataset = SyntheticDataset(batch_size=100, nb_genes=100, n_batches=3)

        n_genes = dataset.nb_genes
        n_top = n_genes // 2
        dataset._highly_variable_genes(n_bins=3, flavor="seurat_v2")
        df = dataset._highly_variable_genes(
            n_bins=3, n_top_genes=n_top, flavor="seurat_v2"
        )
        assert df["highly_variable"].sum() >= n_top

        dataset.filter_genes_by_count(2, per_batch=True)
        dataset.subsample_genes(new_n_genes=n_top)
        new_genes = dataset.nb_genes
        assert n_genes > new_genes, "subsample_genes did not filter out genes"

        dataset = SyntheticDataset(batch_size=100, nb_genes=100, n_batches=3)
        n_genes = dataset.nb_genes
        n_top = n_genes // 2
        df = dataset._highly_variable_genes(
            n_bins=3, flavor="seurat_v2", batch_correction=False, n_top_genes=n_top
        )
        assert (
            "highly_variable_nbatches" not in df.columns
        ), "HVG dataframe should not contain batch information"
        df = dataset._highly_variable_genes(
            n_bins=3, flavor="seurat_v2", batch_correction=True
        )
        assert "highly_variable_nbatches" in df.columns
        assert "highly_variable_intersection" in df.columns
        df = dataset._highly_variable_genes(
            n_bins=3, flavor="seurat_v3", batch_correction=False, n_top_genes=n_top
        )
        assert (
            "highly_variable_nbatches" not in df.columns
        ), "HVG dataframe should not contain batch information"
        df = dataset._highly_variable_genes(
            n_bins=3, flavor="seurat_v3", batch_correction=True, n_top_genes=n_top
        )
        assert "highly_variable_nbatches" in df.columns
        assert "highly_variable_intersection" in df.columns

    def test_dense_subsample_genes(self):
        dataset = SyntheticDataset(batch_size=100, nb_genes=100, n_batches=3)

        n_genes = dataset.nb_genes
        n_top = n_genes // 2
        dataset.subsample_genes(new_n_genes=n_top, mode="cell_ranger")
        assert dataset.nb_genes == n_top

        # With Seurat v2
        dataset = SyntheticDataset(batch_size=100, nb_genes=100, n_batches=3)
        dataset.subsample_genes(new_n_genes=n_top, mode="seurat_v2")
        assert dataset.nb_genes == n_top

        # With Seurat v3
        dataset = SyntheticDataset(batch_size=100, nb_genes=100, n_batches=3)
        dataset.subsample_genes(new_n_genes=n_top, mode="seurat_v3")
        assert dataset.nb_genes == n_top

        # make sure constant genes have low scores
        dataset = SyntheticDataset(batch_size=100, nb_genes=100, n_batches=3)
        dataset.X[:, -1] = np.zeros_like(dataset.X[:, -1])
        df = dataset._highly_variable_genes(n_top_genes=n_top, flavor="seurat_v3")

        assert df.loc[str(dataset.nb_genes - 1)]["highly_variable_median_variance"] == 0
