import pytest


def test_poisson_gene_selection():
    import numpy as np
    import pytest

    from scvi.data import poisson_gene_selection, synthetic_iid

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
    adata = synthetic_iid(sparse_format="csr_matrix")
    poisson_gene_selection(adata, batch_key="batch", n_top_genes=n_top_genes)
    assert np.sum(adata.var["highly_variable"]) == n_top_genes

    X = adata.X
    adata.X = -X
    with pytest.raises(ValueError):
        poisson_gene_selection(adata, batch_key="batch", n_top_genes=n_top_genes)
    adata.X = 0.25 * X
    with pytest.raises(ValueError):
        poisson_gene_selection(adata, batch_key="batch", n_top_genes=n_top_genes)


@pytest.mark.internet
def test_add_dna_sequence(save_path: str):
    from scvi.data import add_dna_sequence, synthetic_iid

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
        genome_dir=save_path,
        install_genome=False,
    )
    assert "dna_sequence" in adata.varm.keys()
    assert "dna_code" in adata.varm.keys()
    assert adata.varm["dna_code"].values.shape[1] == seq_len


def test_reads_to_fragments():
    from scvi.data import reads_to_fragments, synthetic_iid

    adata = synthetic_iid()
    reads_to_fragments(adata)

    assert "fragments" in adata.layers.keys()

    adata = synthetic_iid(sparse_format="csr_matrix")
    reads_to_fragments(adata)

    assert "fragments" in adata.layers.keys()

    adata = synthetic_iid()
    adata.layers["reads"] = adata.X
    reads_to_fragments(adata, read_layer="reads")

    assert "fragments" in adata.layers.keys()
