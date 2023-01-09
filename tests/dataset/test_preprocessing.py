import numpy as np
import pytest
from scipy import sparse

from scvi.data import add_dna_sequence, poisson_gene_selection, synthetic_iid


def test_poisson_gene_selection():
    n_top_genes = 10
    adata = synthetic_iid()
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


@pytest.mark.internet
def test_add_dna_sequence(save_path):
    adata = synthetic_iid()
    adata = adata[:, :2].copy()
    adata.var["chr"] = "chr1"
    adata.var["start"] = [629395, 633578]
    adata.var["start"] = [630394, 634591]

    add_dna_sequence(adata, seq_len=100, genome="hg38", genome_dir=save_path)
    assert "dna_sequence" in adata.varm.keys()
    assert "dna_code" in adata.varm.keys()
    assert adata.varm["dna_code"].values.shape[1] == 100
