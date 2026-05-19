import tempfile

import numpy as np
import pandas as pd
import pytest
from anndata import AnnData

import scvi.external.harreman.hotspot as hs
import scvi.external.harreman.preprocessing as pp
import scvi.external.harreman.tools as tl
from scvi.external.harreman.preprocessing.anndata import counts_from_anndata


@pytest.fixture
def adata_spatial():
    n_obs = 50
    n_vars = 20
    np.random.seed(42)
    X = np.random.poisson(1.0, size=(n_obs, n_vars)).astype(float)
    obs = pd.DataFrame(
        {"sample": pd.Categorical(np.random.choice(["s1", "s2"], size=n_obs))},
        index=[f"cell{i}" for i in range(n_obs)],
    )
    var = pd.DataFrame(index=[f"gene{i}" for i in range(n_vars)])
    obsm = {"spatial": np.random.rand(n_obs, 2) * 100}
    return AnnData(X=X, obs=obs, var=var, obsm=obsm)


@pytest.fixture
def adata_with_modules(adata_spatial):
    """AnnData with KNN graph, autocorrelation and local correlation computed."""
    tl.compute_knn_graph(adata_spatial, compute_neighbors_on_key="spatial", n_neighbors=5)
    hs.compute_local_autocorrelation(adata_spatial, model="danb")
    genes = adata_spatial.var_names[:10].tolist()
    hs.compute_local_correlation(adata_spatial, genes=genes)
    return adata_spatial


@pytest.fixture
def adata_with_graph(adata_spatial):
    tl.compute_knn_graph(adata_spatial, compute_neighbors_on_key="spatial", n_neighbors=5)
    return adata_spatial


@pytest.fixture
def adata_with_autocorrelation(adata_with_graph):
    hs.compute_local_autocorrelation(adata_with_graph, model="danb")
    return adata_with_graph


# ── KNN graph ─────────────────────────────────────────────────────────────────


def test_compute_knn_graph_n_neighbors(adata_spatial):
    tl.compute_knn_graph(adata_spatial, compute_neighbors_on_key="spatial", n_neighbors=5)
    assert "distances" in adata_spatial.obsp
    assert "weights" in adata_spatial.obsp
    assert adata_spatial.obsp["distances"].shape == (adata_spatial.n_obs, adata_spatial.n_obs)


def test_compute_knn_graph_weighted(adata_spatial):
    tl.compute_knn_graph(
        adata_spatial, compute_neighbors_on_key="spatial", n_neighbors=5, weighted_graph=True
    )
    assert "weights" in adata_spatial.obsp
    assert adata_spatial.obsp["weights"].nnz > 0


def test_compute_knn_graph_neighborhood_radius(adata_spatial):
    tl.compute_knn_graph(
        adata_spatial, compute_neighbors_on_key="spatial", neighborhood_radius=20.0
    )
    assert "distances" in adata_spatial.obsp
    assert "weights" in adata_spatial.obsp


def test_compute_knn_graph_with_sample_key(adata_spatial):
    tl.compute_knn_graph(
        adata_spatial, compute_neighbors_on_key="spatial", n_neighbors=5, sample_key="sample"
    )
    assert "distances" in adata_spatial.obsp
    assert adata_spatial.uns["sample_key"] == "sample"


def test_compute_knn_graph_from_distances(adata_spatial):
    tl.compute_knn_graph(adata_spatial, compute_neighbors_on_key="spatial", n_neighbors=5)
    adata2 = adata_spatial.copy()
    del adata2.obsp["weights"]
    tl.compute_knn_graph(adata2, distances_obsp_key="distances", n_neighbors=5)
    assert "weights" in adata2.obsp


def test_compute_knn_graph_missing_key(adata_spatial):
    with pytest.raises(ValueError, match="not found in adata.obsm"):
        tl.compute_knn_graph(
            adata_spatial, compute_neighbors_on_key="nonexistent_key", n_neighbors=5
        )


def test_compute_knn_graph_no_neighbors_raises(adata_spatial):
    with pytest.raises(ValueError, match="Either 'n_neighbors' or 'neighborhood_radius'"):
        tl.compute_knn_graph(adata_spatial, compute_neighbors_on_key="spatial")


# ── Local autocorrelation ──────────────────────────────────────────────────────


def test_compute_local_autocorrelation(adata_with_graph):
    hs.compute_local_autocorrelation(adata_with_graph, model="danb")
    assert "gene_autocorrelation_results" in adata_with_graph.uns
    results = adata_with_graph.uns["gene_autocorrelation_results"]
    assert len(results) == adata_with_graph.n_vars


@pytest.mark.parametrize("model", ["none", "normal", "bernoulli"])
def test_compute_local_autocorrelation_models(adata_with_graph, model):
    hs.compute_local_autocorrelation(adata_with_graph, model=model)
    results = adata_with_graph.uns["gene_autocorrelation_results"]
    assert len(results) == adata_with_graph.n_vars


def test_compute_local_autocorrelation_result_columns(adata_with_graph):
    hs.compute_local_autocorrelation(adata_with_graph, model="danb")
    cols = adata_with_graph.uns["gene_autocorrelation_results"].columns.tolist()
    assert cols == ["C", "Z", "Z_Pval", "Z_FDR"]


def test_compute_local_autocorrelation_with_explicit_genes(adata_with_graph):
    genes = adata_with_graph.var_names[:8].tolist()
    hs.compute_local_autocorrelation(adata_with_graph, model="danb", genes=genes)
    results = adata_with_graph.uns["gene_autocorrelation_results"]
    assert len(results) == 8
    assert set(results.index).issubset(set(genes))


def test_compute_local_autocorrelation_sample_specific(adata_spatial):
    tl.compute_knn_graph(
        adata_spatial, compute_neighbors_on_key="spatial", n_neighbors=5, sample_key="sample"
    )
    hs.compute_local_autocorrelation(adata_spatial, model="danb")
    assert "gene_autocorrelation_results" in adata_spatial.uns
    assert len(adata_spatial.uns["gene_autocorrelation_results"]) == adata_spatial.n_vars


def test_compute_local_autocorrelation_permutation_test(adata_with_graph):
    hs.compute_local_autocorrelation(adata_with_graph, model="danb", permutation_test=True, M=5)
    cols = adata_with_graph.uns["gene_autocorrelation_results"].columns.tolist()
    assert "Perm_Pval" in cols
    assert "Perm_FDR" in cols


def test_compute_local_autocorrelation_sorted_by_z(adata_with_graph):
    hs.compute_local_autocorrelation(adata_with_graph, model="danb")
    z_vals = adata_with_graph.uns["gene_autocorrelation_results"]["Z"].values
    assert (z_vals[:-1] >= z_vals[1:]).all()


# ── Local correlation ──────────────────────────────────────────────────────────


def test_compute_local_correlation(adata_with_autocorrelation):
    genes = adata_with_autocorrelation.var_names[:5].tolist()
    hs.compute_local_correlation(adata_with_autocorrelation, genes=genes)
    assert "lc_zs" in adata_with_autocorrelation.uns


def test_compute_local_correlation_result_keys(adata_with_autocorrelation):
    genes = adata_with_autocorrelation.var_names[:5].tolist()
    hs.compute_local_correlation(adata_with_autocorrelation, genes=genes)
    for key in ["lcs", "lc_zs", "lc_z_pvals", "lc_z_FDR"]:
        assert key in adata_with_autocorrelation.uns


def test_compute_local_correlation_result_shape(adata_with_autocorrelation):
    genes = adata_with_autocorrelation.var_names[:5].tolist()
    hs.compute_local_correlation(adata_with_autocorrelation, genes=genes)
    lc_zs = adata_with_autocorrelation.uns["lc_zs"]
    assert lc_zs.shape == (5, 5)
    assert list(lc_zs.index) == genes
    assert list(lc_zs.columns) == genes


def test_compute_local_correlation_permutation_test(adata_with_autocorrelation):
    genes = adata_with_autocorrelation.var_names[:4].tolist()
    hs.compute_local_correlation(
        adata_with_autocorrelation, genes=genes, permutation_test=True, M=5
    )
    assert "lc_perm_pvals" in adata_with_autocorrelation.uns
    assert "lc_perm_pvals_sym" in adata_with_autocorrelation.uns


# ── Modules ────────────────────────────────────────────────────────────────────


def test_create_modules(adata_with_autocorrelation):
    genes = adata_with_autocorrelation.var_names.tolist()
    hs.compute_local_correlation(adata_with_autocorrelation, genes=genes)
    hs.create_modules(adata_with_autocorrelation, min_gene_threshold=2)
    for key in ["modules", "gene_modules_dict", "linkage"]:
        assert key in adata_with_autocorrelation.uns


def test_calculate_module_scores(adata_with_autocorrelation):
    genes = adata_with_autocorrelation.var_names.tolist()
    hs.compute_local_correlation(adata_with_autocorrelation, genes=genes)
    hs.create_modules(adata_with_autocorrelation, min_gene_threshold=2)
    modules = adata_with_autocorrelation.uns["gene_modules_dict"]
    if not any(k != "-1" for k in modules):
        pytest.skip("No modules found in random data")
    hs.calculate_module_scores(adata_with_autocorrelation)
    assert "module_scores" in adata_with_autocorrelation.obsm
    assert "gene_loadings" in adata_with_autocorrelation.varm


def test_compute_top_scoring_modules(adata_with_autocorrelation):
    genes = adata_with_autocorrelation.var_names.tolist()
    hs.compute_local_correlation(adata_with_autocorrelation, genes=genes)
    hs.create_modules(adata_with_autocorrelation, min_gene_threshold=2)
    modules = adata_with_autocorrelation.uns["gene_modules_dict"]
    if not any(k != "-1" for k in modules):
        pytest.skip("No modules found in random data")
    hs.calculate_module_scores(adata_with_autocorrelation)
    top = hs.compute_top_scoring_modules(adata_with_autocorrelation)
    assert len(top) == adata_with_autocorrelation.n_obs


# ── Preprocessing ──────────────────────────────────────────────────────────────


def test_counts_from_anndata_default_layer(adata_spatial):
    counts = counts_from_anndata(adata_spatial, dense=True)
    assert counts.shape == (adata_spatial.n_vars, adata_spatial.n_obs)


def test_counts_from_anndata_explicit_layer(adata_spatial):
    adata_spatial.layers["counts"] = adata_spatial.X.copy()
    counts = counts_from_anndata(adata_spatial, layer_key="counts", dense=True)
    assert counts.shape == (adata_spatial.n_vars, adata_spatial.n_obs)


def test_write_read_h5ad_roundtrip(adata_with_autocorrelation):
    genes = adata_with_autocorrelation.var_names[:5].tolist()
    hs.compute_local_correlation(adata_with_autocorrelation, genes=genes)
    with tempfile.NamedTemporaryFile(suffix=".h5ad", delete=False) as f:
        path = f.name
    pp.write_h5ad(adata_with_autocorrelation, filename=path)
    recovered = pp.read_h5ad(path)
    assert "gene_autocorrelation_results" in recovered.uns
    assert "lc_zs" in recovered.uns


def test_write_h5ad_no_filename_raises(adata_spatial):
    with pytest.raises(ValueError, match="path to save"):
        pp.write_h5ad(adata_spatial)
