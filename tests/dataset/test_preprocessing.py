import os

import numpy as np
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


def test_add_dna_sequence():
    adata = synthetic_iid()
    adata = adata[:, :2].copy()
    adata.var["chr"] = "chr1"
    adata.var["start"] = [2, 20]
    adata.var["end"] = [30, 80]
    seq_len = 6
    add_dna_sequence(
        adata,
        seq_len=seq_len,
        genome_name="test_genome",
        genome_dir=os.path.abspath("tests/data/"),
        install_genome=False,
    )
    assert "dna_sequence" in adata.varm.keys()
    assert "dna_code" in adata.varm.keys()
    assert adata.varm["dna_code"].values.shape[1] == seq_len
