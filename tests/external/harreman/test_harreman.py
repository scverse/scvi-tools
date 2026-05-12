import numpy as np
import pandas as pd
import pytest
from anndata import AnnData

import scvi.external.harreman.tools as tl
import scvi.external.harreman.hotspot as hs


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


def test_compute_knn_graph_n_neighbors(adata_spatial):
    tl.compute_knn_graph(adata_spatial, compute_neighbors_on_key="spatial", n_neighbors=5)
    assert "distances" in adata_spatial.obsp
    assert "weights" in adata_spatial.obsp
    assert adata_spatial.obsp["distances"].shape == (adata_spatial.n_obs, adata_spatial.n_obs)


def test_compute_knn_graph_weighted(adata_spatial):
    tl.compute_knn_graph(adata_spatial, compute_neighbors_on_key="spatial", n_neighbors=5, weighted_graph=True)
    assert "weights" in adata_spatial.obsp
    assert adata_spatial.obsp["weights"].nnz > 0


def test_compute_knn_graph_with_sample_key(adata_spatial):
    tl.compute_knn_graph(adata_spatial, compute_neighbors_on_key="spatial", n_neighbors=5, sample_key="sample")
    assert "distances" in adata_spatial.obsp


def test_compute_knn_graph_missing_key(adata_spatial):
    with pytest.raises(ValueError, match="not found in adata.obsm"):
        tl.compute_knn_graph(adata_spatial, compute_neighbors_on_key="nonexistent_key", n_neighbors=5)


def test_compute_knn_graph_no_neighbors_raises(adata_spatial):
    with pytest.raises(ValueError, match="Either \'n_neighbors\' or \'neighborhood_radius\'"):
        tl.compute_knn_graph(adata_spatial, compute_neighbors_on_key="spatial")


def test_compute_local_autocorrelation(adata_spatial):
    tl.compute_knn_graph(adata_spatial, compute_neighbors_on_key="spatial", n_neighbors=5)
    hs.compute_local_autocorrelation(adata_spatial, model="danb")
    assert "gene_autocorrelation_results" in adata_spatial.uns
    results = adata_spatial.uns["gene_autocorrelation_results"]
    assert len(results) == adata_spatial.n_vars


def test_compute_local_correlation(adata_spatial):
    tl.compute_knn_graph(adata_spatial, compute_neighbors_on_key="spatial", n_neighbors=5)
    hs.compute_local_autocorrelation(adata_spatial, model="danb")
    # Pass genes explicitly to avoid empty selection from FDR filtering
    genes = adata_spatial.var_names[:5].tolist()
    hs.compute_local_correlation(adata_spatial, genes=genes)
    assert "lc_zs" in adata_spatial.uns
