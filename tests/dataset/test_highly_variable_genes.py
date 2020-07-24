import unittest
import numpy as np
from unittest import TestCase
from scvi.dataset import synthetic_iid
from scvi.dataset import highly_variable_genes_seurat_v3  # , poisson_gene_selection


class TestHighlyVariableGenes(TestCase):
    def test_highly_variable_genes_seurat_v3(self):
        adata = synthetic_iid()

        n_top_genes = 20
        highly_variable_genes_seurat_v3(adata, n_top_genes=n_top_genes)
        keys = [
            "highly_variable",
            "highly_variable_rank",
            "means",
            "variances",
            "variances_norm",
        ]
        for key in keys:
            assert key in adata.var.keys()
        assert np.sum(adata.var["highly_variable"]) == n_top_genes
        assert np.max(adata.var["highly_variable_rank"]) == n_top_genes - 1
        highly_variable_genes_seurat_v3(
            adata, n_top_genes=n_top_genes, batch_key="batch"
        )

    def test_poisson_gene_selection(self):
        pass


if __name__ == "__main__":
    unittest.main()
