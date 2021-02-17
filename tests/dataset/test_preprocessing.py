import numpy as np
from scipy import sparse

from scvi.data import poisson_gene_selection, synthetic_iid


def test_poisson_gene_selection():
    adata = synthetic_iid()
    n_top_genes = 10
    poisson_gene_selection(adata, batch_key="batch", n_top_genes=n_top_genes)
    keys = [
        "highly_variable",
        "observed_fraction_zeros",
        "expected_fraction_zeros",
        "prob_zero_enriched_nbatches",
        "prob_zero_enrichment",
        "prob_zero_enrichment_rank",
    ]
    for key in keys:
        assert key in adata.var.keys()
    assert np.sum(adata.var["highly_variable"]) == n_top_genes
    adata = synthetic_iid()
    adata.X = sparse.csr_matrix(adata.X)
    poisson_gene_selection(adata, batch_key="batch", n_top_genes=n_top_genes)
    assert np.sum(adata.var["highly_variable"]) == n_top_genes
