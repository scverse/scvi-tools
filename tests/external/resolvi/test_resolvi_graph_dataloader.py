"""GraphDataLoader integration tests for RESOLVI."""

import time
import warnings

import numpy as np
import pytest
import torch

from scvi.data import synthetic_iid
from scvi.external import RESOLVI


@pytest.fixture
def adata():
    adata = synthetic_iid(generate_coordinates=True, n_regions=5, n_proteins=10)
    adata.obsm["X_spatial"] = adata.obsm["coordinates"]
    adata.obs["cell_area"] = np.random.default_rng(0).gamma(2.0, 1.0, size=adata.n_obs)
    RESOLVI.setup_anndata(adata)
    return adata


def _resolvi_graph_cls():
    from scvi.dataloaders import GraphDataSplitter

    class RESOLVIGraph(RESOLVI):
        _data_splitter_cls = GraphDataSplitter

    return RESOLVIGraph


def _resolvi_legacy_cls():
    from scvi.dataloaders import DataSplitter

    class RESOLVILegacy(RESOLVI):
        _data_splitter_cls = DataSplitter

    return RESOLVILegacy


def _train_graph(model, max_epochs: int = 2, **kwargs):
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        model.train(
            max_epochs=max_epochs,
            datasplitter_kwargs={"neighbor_indices_key": "index_neighbor"},
            **kwargs,
        )


def test_resolvi_get_fn_args_prefers_graph_batch_x_n(adata):
    """RESOLVI must use pre-fetched GraphDataLoader neighbor expression when present."""
    from scvi.dataloaders import GraphDataLoader

    adata_manager = RESOLVI._get_most_recent_anndata_manager(adata, required=True)
    model = RESOLVI(adata)
    batch = next(
        iter(
            GraphDataLoader(
                adata_manager,
                full_adata_manager=adata_manager,
                batch_size=8,
                shuffle=False,
            )
        )
    )

    class FailingExpressionDataset:
        def __getitem__(self, item):
            raise AssertionError("fallback AnnTorchDataset should not be used")

    model.module.model.expression_anntorchdata = FailingExpressionDataset()
    _, kwargs = model.module._get_fn_args_from_batch(batch)

    torch.testing.assert_close(kwargs["x_n"], batch.x_n.reshape(batch.x.shape[0], -1))
    torch.testing.assert_close(kwargs["distances_n"], batch.distances_n)


def test_resolvi_graph_path_trains(adata):
    """RESOLVI must train to completion when forced to use GraphDataLoader."""
    RESOLVIGraph = _resolvi_graph_cls()

    RESOLVIGraph.setup_anndata(adata)
    model = RESOLVIGraph(adata)
    _train_graph(model)

    latent = model.get_latent_representation()
    assert latent.shape == (adata.n_obs, model.module.n_latent)

    model = RESOLVIGraph(adata, dispersion="gene-batch")
    _train_graph(model)


def test_resolvi_uses_graph_datasplitter_by_default():
    """RESOLVI should opt into GraphDataSplitter without a test-only subclass."""
    from scvi.dataloaders import GraphDataSplitter

    assert RESOLVI._data_splitter_cls is GraphDataSplitter


def test_resolvi_graph_train_size_factor(adata):
    """GraphDataLoader path supports RESOLVI size-factor training modes."""
    RESOLVIGraph = _resolvi_graph_cls()

    RESOLVIGraph.setup_anndata(adata, batch_key="batch", size_factor_key="cell_area")
    model = RESOLVIGraph(adata, size_scaling=True)
    _train_graph(model)

    RESOLVIGraph.setup_anndata(adata, size_factor_key="cell_area")
    model = RESOLVIGraph(adata, size_scaling=False)
    _train_graph(model)


@pytest.mark.optional
def test_resolvi_graph_save_load(adata, tmp_path):
    """GraphDataLoader path preserves legacy save/load behavior."""
    RESOLVIGraph = _resolvi_graph_cls()

    RESOLVIGraph.setup_anndata(adata)
    model = RESOLVIGraph(adata)
    _train_graph(model)
    hist_elbo = model.history_["elbo_train"]
    latent = model.get_latent_representation()
    assert latent.shape == (adata.n_obs, model.module.n_latent)
    model.differential_expression(groupby="labels")
    model.differential_expression(groupby="labels", weights="importance")

    save_path = str(tmp_path / "test_resolvi_graph")
    model.save(save_path, save_anndata=True, overwrite=True)
    model2 = model.load(save_path)
    np.testing.assert_array_equal(model2.history_["elbo_train"], hist_elbo)
    latent2 = model2.get_latent_representation()
    assert np.allclose(latent, latent2)
    model.load_query_data(reference_model=save_path, adata=adata)


@pytest.mark.optional
def test_resolvi_graph_downstream(adata, tmp_path):
    """GraphDataLoader path covers legacy downstream RESOLVI APIs."""
    RESOLVIGraph = _resolvi_graph_cls()

    RESOLVIGraph.setup_anndata(adata, size_factor_key="cell_area")
    model = RESOLVIGraph(adata)
    _train_graph(model)
    latent = model.get_latent_representation()
    assert latent.shape == (adata.n_obs, model.module.n_latent)
    _ = model.get_normalized_expression(n_samples=31, library_size=10000)
    _ = model.get_normalized_expression_importance(n_samples=30, library_size=10000)
    _ = model.get_normalized_expression_importance(n_samples=30, size_scaling=True)
    model.differential_expression(groupby="labels")
    model.differential_expression(groupby="labels", weights="importance")
    model.differential_expression(groupby="labels", weights="importance", size_scaling=True)
    model.sample_posterior(
        model=model.module.model_residuals,
        num_samples=30,
        return_samples=False,
        return_sites=None,
        batch_size=1000,
    )
    model.sample_posterior(
        model=model.module.model_residuals,
        num_samples=30,
        return_samples=False,
        batch_size=1000,
    )

    model.load_query_data(reference_model=model, adata=adata)
    save_path = str(tmp_path / "test_resolvi_graph")
    model.save(save_path, save_anndata=True, overwrite=True)
    model_query = model.load_query_data(reference_model=save_path, adata=adata)
    _train_graph(model_query)


def test_resolvi_graph_downstream_size_scaling(adata, tmp_path):
    """GraphDataLoader path covers downstream APIs with size scaling enabled."""
    RESOLVIGraph = _resolvi_graph_cls()

    RESOLVIGraph.setup_anndata(adata, size_factor_key="cell_area")
    model = RESOLVIGraph(adata, size_scaling=True)
    _train_graph(model)
    latent = model.get_latent_representation()
    assert latent.shape == (adata.n_obs, model.module.n_latent)
    _ = model.get_normalized_expression(n_samples=31, library_size=10000)
    _ = model.get_normalized_expression_importance(n_samples=30, library_size=10000)
    _ = model.get_normalized_expression_importance(n_samples=30, size_scaling=True)
    model.differential_expression(groupby="labels")
    model.differential_expression(groupby="labels", weights="importance")
    model.differential_expression(groupby="labels", weights="importance", size_scaling=True)
    model.sample_posterior(
        model=model.module.model_residuals,
        num_samples=30,
        return_samples=False,
        return_sites=None,
        batch_size=1000,
    )
    model.sample_posterior(
        model=model.module.model_residuals,
        num_samples=30,
        return_samples=False,
        batch_size=1000,
    )

    model.load_query_data(reference_model=model, adata=adata)
    save_path = str(tmp_path / "test_resolvi_graph")
    model.save(save_path, save_anndata=True, overwrite=True)
    model_query = model.load_query_data(reference_model=save_path, adata=adata)
    _train_graph(model_query)


@pytest.mark.optional
def test_resolvi_graph_semisupervised(adata):
    """GraphDataLoader path supports semisupervised RESOLVI APIs."""
    RESOLVIGraph = _resolvi_graph_cls()

    RESOLVIGraph.setup_anndata(adata, labels_key="labels")
    model = RESOLVIGraph(adata, semisupervised=True)
    _train_graph(model)
    model.differential_niche_abundance(
        batch_size=30,
        groupby="batch",
        neighbor_key="index_neighbor",
    )
    pred = model.predict(soft=True)
    assert pred.shape == (adata.n_obs, model.summary_stats.n_labels - 1)
    pred = model.predict(soft=False)
    assert pred.shape == (adata.n_obs,)


def test_resolvi_graph_scarches(adata):
    """GraphDataLoader path preserves legacy scArches query workflow."""
    RESOLVIGraph = _resolvi_graph_cls()

    adata.obs["hemisphere"] = ["right" if x > 0 else "left" for x in adata.obsm["X_spatial"][:, 0]]
    ref_adata = adata[adata.obs["hemisphere"] == "left"].copy()
    query_adata = adata[adata.obs["hemisphere"] == "right"].copy()

    RESOLVIGraph.setup_anndata(ref_adata, labels_key="labels")
    model = RESOLVIGraph(ref_adata, semisupervised=True)
    _train_graph(model)

    ref_adata.obsm["resolvi_celltypes"] = model.predict(ref_adata, num_samples=3, soft=True)
    ref_adata.obs["resolvi_predicted"] = ref_adata.obsm["resolvi_celltypes"].idxmax(axis=1)
    ref_adata.obsm["X_resolVI"] = model.get_latent_representation(ref_adata)

    query_adata.obs["predicted_celltype"] = "unknown"
    query_adata.obs_names = [f"query_{i}" for i in query_adata.obs_names]

    model.prepare_query_anndata(query_adata, reference_model=model)
    query_resolvi = model.load_query_data(query_adata, reference_model=model)

    _train_graph(query_resolvi, max_epochs=1)

    query_adata.obs["resolvi_predicted"] = query_resolvi.predict(
        query_adata,
        num_samples=3,
        soft=False,
    )
    query_adata.obsm["X_resolVI"] = query_resolvi.get_latent_representation(query_adata)


@pytest.mark.parametrize("weights", ["importance", "uniform"])
@pytest.mark.parametrize("n_samples", [1, 3])
@pytest.mark.parametrize("downsample_counts", [True, False])
def test_resolvi_graph_differential_expression(
    adata,
    weights: str,
    n_samples: int,
    downsample_counts: bool,
):
    """GraphDataLoader path supports RESOLVI differential expression settings."""
    RESOLVIGraph = _resolvi_graph_cls()

    RESOLVIGraph.setup_anndata(adata)
    model = RESOLVIGraph(adata, downsample_counts=downsample_counts)
    _train_graph(model, max_epochs=1)
    model.differential_expression(groupby="labels", weights=weights, n_samples=n_samples)


@pytest.mark.benchmark
def test_resolvi_dataloader_speed_comparison(adata):
    """Side-by-side wall-clock comparison: AnnDataLoader vs GraphDataLoader."""
    RESOLVILegacy = _resolvi_legacy_cls()
    RESOLVIGraph = _resolvi_graph_cls()

    n_epochs = 5

    RESOLVILegacy.setup_anndata(adata)
    model_ann = RESOLVILegacy(adata)
    t0 = time.perf_counter()
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        model_ann.train(max_epochs=n_epochs)
    t_ann = time.perf_counter() - t0

    RESOLVIGraph.setup_anndata(adata)
    model_graph = RESOLVIGraph(adata)
    t0 = time.perf_counter()
    _train_graph(model_graph, max_epochs=n_epochs)
    t_graph = time.perf_counter() - t0

    print(f"\nAnnDataLoader:   {t_ann:.2f}s total  ({t_ann / n_epochs:.3f}s/epoch)")
    print(f"GraphDataLoader: {t_graph:.2f}s total  ({t_graph / n_epochs:.3f}s/epoch)")
    print(f"Ratio (graph/ann): {t_graph / t_ann:.2f}x")

    assert t_graph / t_ann < 3.0, (
        f"GraphDataLoader is {t_graph / t_ann:.1f}x slower than AnnDataLoader"
    )


@pytest.mark.benchmark
def test_resolvi_elbo_comparable_between_paths(adata):
    """Both paths should reach similar final ELBO after training."""
    RESOLVILegacy = _resolvi_legacy_cls()
    RESOLVIGraph = _resolvi_graph_cls()

    n_epochs = 10

    RESOLVILegacy.setup_anndata(adata)
    model_ann = RESOLVILegacy(adata)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        model_ann.train(max_epochs=n_epochs)
    elbo_ann = model_ann.history_["elbo_train"].iloc[-1].values[0]

    RESOLVIGraph.setup_anndata(adata)
    model_graph = RESOLVIGraph(adata)
    _train_graph(model_graph, max_epochs=n_epochs)
    elbo_graph = model_graph.history_["elbo_train"].iloc[-1].values[0]

    print(f"\nFinal ELBO - AnnDataLoader: {elbo_ann:.2f} GraphDataLoader: {elbo_graph:.2f}")

    assert np.isfinite(elbo_graph), "GraphDataLoader ELBO is not finite"
    assert abs(elbo_graph - elbo_ann) / (abs(elbo_ann) + 1e-8) < 0.5, (
        f"ELBO diverged: ann={elbo_ann:.2f} graph={elbo_graph:.2f}"
    )


@pytest.mark.benchmark
def test_resolvi_graph_elbo_decreases(adata):
    """ELBO must decrease over training with GraphDataLoader."""
    RESOLVIGraph = _resolvi_graph_cls()

    RESOLVIGraph.setup_anndata(adata)
    model = RESOLVIGraph(adata)
    _train_graph(model, max_epochs=10)

    history = model.history_["elbo_train"]
    assert history.iloc[-1].values[0] < history.iloc[0].values[0], (
        "ELBO did not decrease with GraphDataLoader"
    )


@pytest.mark.benchmark
def test_resolvi_latent_shape_graph_path(adata):
    """get_latent_representation() must return (n_obs, n_latent) with GraphDataLoader."""
    RESOLVIGraph = _resolvi_graph_cls()

    RESOLVIGraph.setup_anndata(adata)
    model = RESOLVIGraph(adata)
    _train_graph(model, max_epochs=3)
    latent = model.get_latent_representation()
    assert latent.shape == (adata.n_obs, model.module.n_latent)
