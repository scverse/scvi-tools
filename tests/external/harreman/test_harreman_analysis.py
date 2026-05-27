from __future__ import annotations

import os
import tempfile
from unittest.mock import MagicMock, patch

import numba
import numpy as np
import pandas as pd
import pytest
from anndata import AnnData

from scvi.external.harreman._analysis import HarremanAnalysis
from scvi.external.harreman._constants import (
    HARREMAN_AUTOCORR_KEY,
    HARREMAN_CCC_KEY,
    HARREMAN_CT_CCC_KEY,
    HARREMAN_GENE_PAIRS_KEY,
    HARREMAN_ICS_KEY,
    HARREMAN_PARAMS_KEY,
    HARREMAN_SIG_KEY,
    HARREMAN_UNS_KEY,
    STEP_CCC,
    STEP_FILTER,
    STEP_GENE_PAIRS,
    STEP_ICS,
    STEP_SETUP,
    STEP_SIG,
)
from scvi.external.harreman._results import HarremanResults

_NUMBA_CACHE_DIR = os.path.join(tempfile.gettempdir(), "numba_cache")
os.environ.setdefault("NUMBA_CACHE_DIR", _NUMBA_CACHE_DIR)
numba.config.CACHE_DIR = _NUMBA_CACHE_DIR

# ── Shared fixtures ────────────────────────────────────────────────────────────


@pytest.fixture
def adata_spatial():
    n_obs, n_vars = 50, 20
    rng = np.random.default_rng(42)
    X = rng.poisson(1.0, size=(n_obs, n_vars)).astype(float)
    obs = pd.DataFrame(
        {"cell_type": pd.Categorical(np.tile(["TypeA", "TypeB"], n_obs // 2))},
        index=[f"cell{i}" for i in range(n_obs)],
    )
    var = pd.DataFrame(index=[f"gene{i}" for i in range(n_vars)])
    obsm = {"spatial": rng.random((n_obs, 2)) * 100}
    return AnnData(X=X, obs=obs, var=var, obsm=obsm)


@pytest.fixture
def uns_with_results():
    df = pd.DataFrame({"score": [0.1, 0.2]})
    return {
        HARREMAN_AUTOCORR_KEY: df.copy(),
        HARREMAN_GENE_PAIRS_KEY: df.copy(),
        HARREMAN_CCC_KEY: df.copy(),
        HARREMAN_CT_CCC_KEY: None,
        HARREMAN_ICS_KEY: None,
        HARREMAN_SIG_KEY: df.copy(),
        HARREMAN_PARAMS_KEY: {"species": "human"},
    }


# ── HarremanResults tests ──────────────────────────────────────────────────────


def test_results_from_uns(uns_with_results):
    r = HarremanResults.from_uns(uns_with_results)
    assert isinstance(r.autocorrelation, pd.DataFrame)
    assert isinstance(r.cell_communication, pd.DataFrame)
    assert r.ct_cell_communication is None
    assert r.params == {"species": "human"}


def test_results_from_uns_missing_params_key_raises(uns_with_results):
    del uns_with_results[HARREMAN_PARAMS_KEY]
    with pytest.raises(KeyError):
        HarremanResults.from_uns(uns_with_results)


def test_results_repr_does_not_raise(uns_with_results):
    r = HarremanResults.from_uns(uns_with_results)
    repr(r)  # should not raise


# ── Init / basic validation ────────────────────────────────────────────────────


def test_init_stores_adata(adata_spatial):
    ha = HarremanAnalysis(adata_spatial)
    assert ha.adata is adata_spatial


def test_init_non_anndata_raises():
    with pytest.raises(TypeError, match="AnnData"):
        HarremanAnalysis("not_an_adata")


def test_init_is_not_deconvolved_by_default(adata_spatial):
    ha = HarremanAnalysis(adata_spatial)
    assert ha.is_deconvolved is False


def test_init_is_not_set_up_by_default(adata_spatial):
    ha = HarremanAnalysis(adata_spatial)
    assert ha.is_set_up is False


def test_results_before_any_steps_raises(adata_spatial):
    ha = HarremanAnalysis(adata_spatial)
    with pytest.raises(RuntimeError, match="setup"):
        _ = ha.results


def test_filter_genes_before_setup_raises(adata_spatial):
    ha = HarremanAnalysis(adata_spatial)
    with pytest.raises(RuntimeError, match="setup"):
        ha.filter_genes()


def test_compute_gene_pairs_before_setup_raises(adata_spatial):
    ha = HarremanAnalysis(adata_spatial)
    with pytest.raises(RuntimeError, match="setup"):
        ha.compute_gene_pairs()


def test_compute_ccc_before_gene_pairs_raises(adata_spatial):
    ha = HarremanAnalysis(adata_spatial)
    ha._completed_steps.add(STEP_SETUP)  # bypass setup check
    with pytest.raises(RuntimeError, match="compute_gene_pairs"):
        ha.compute_cell_communication()


def test_unsupported_model_raises(adata_spatial):
    bad_model = MagicMock()
    bad_model.__class__.__name__ = "NotAKnownModel"
    with pytest.raises(ValueError, match="Unsupported model"):
        HarremanAnalysis(adata_spatial, model=bad_model)


# ── setup() tests ─────────────────────────────────────────────────────────────


def _make_mock_extract_db(adata):
    """Return a side-effect function that fakes extract_interaction_db without S3."""
    import numpy as np

    def _fake_extract(adata, **kwargs):
        n_vars = adata.n_vars
        db_df = pd.DataFrame(
            np.zeros((n_vars, 3)),
            index=adata.var_names,
            columns=["metab_A", "metab_B", "metab_C"],
        )
        adata.varm["database"] = db_df
        adata.uns["database_varm_key"] = "database"
        adata.uns["database"] = kwargs.get("database", "transporter")
        adata.uns["metabolite_database"] = pd.DataFrame()
        adata.uns["num_metabolites"] = 3
        adata.uns["importer"] = pd.DataFrame()
        adata.uns["exporter"] = pd.DataFrame()
        adata.uns["import_export"] = pd.DataFrame()

    return _fake_extract


def _setup_ha_no_network(adata, monkeypatch, cell_type_key=None):
    """Create HarremanAnalysis with setup() mocked to avoid S3 calls."""
    import scvi.external.harreman._analysis as _mod

    monkeypatch.setattr(_mod, "_extract_interaction_db", _make_mock_extract_db(adata))
    ha = HarremanAnalysis(adata)
    ha.setup(
        compute_neighbors_on_key="spatial",
        n_neighbors=5,
        cell_type_key=cell_type_key,
    )
    return ha


def test_setup_sets_is_set_up(adata_spatial, monkeypatch):
    ha = _setup_ha_no_network(adata_spatial, monkeypatch)
    assert ha.is_set_up is True


def test_setup_builds_knn_graph(adata_spatial, monkeypatch):
    _setup_ha_no_network(adata_spatial, monkeypatch)
    assert "weights" in adata_spatial.obsp


def test_setup_writes_params(adata_spatial, monkeypatch):
    _setup_ha_no_network(adata_spatial, monkeypatch)
    params = adata_spatial.uns[HARREMAN_UNS_KEY][HARREMAN_PARAMS_KEY]
    assert params["species"] == "human"
    assert params["compute_neighbors_on_key"] == "spatial"


def test_setup_with_cell_type_key(adata_spatial, monkeypatch):
    _setup_ha_no_network(adata_spatial, monkeypatch, cell_type_key="cell_type")
    assert adata_spatial.uns["harreman"].get("cell_type_key") == "cell_type"


def test_setup_second_call_overwrites_params(adata_spatial, monkeypatch):
    import scvi.external.harreman._analysis as _mod

    monkeypatch.setattr(_mod, "_extract_interaction_db", _make_mock_extract_db(adata_spatial))
    ha = HarremanAnalysis(adata_spatial)
    ha.setup(compute_neighbors_on_key="spatial", n_neighbors=5, species="human")
    ha.setup(compute_neighbors_on_key="spatial", n_neighbors=5, species="mouse")
    params = adata_spatial.uns[HARREMAN_UNS_KEY][HARREMAN_PARAMS_KEY]
    assert params["species"] == "mouse"


# ── filter_genes() tests ───────────────────────────────────────────────────────


def test_filter_genes_marks_step_complete(adata_spatial, monkeypatch):
    ha = _setup_ha_no_network(adata_spatial, monkeypatch)
    with patch("scvi.external.harreman._analysis._apply_gene_filtering") as mock_filt:
        mock_filt.return_value = None
        ha.filter_genes()
    assert STEP_FILTER in ha._completed_steps


def test_filter_genes_calls_underlying_function(adata_spatial, monkeypatch):
    ha = _setup_ha_no_network(adata_spatial, monkeypatch)
    with patch("scvi.external.harreman._analysis._apply_gene_filtering") as mock_filt:
        mock_filt.return_value = None
        ha.filter_genes(feature_elimination=True, threshold=0.3)
    mock_filt.assert_called_once()
    _, kwargs = mock_filt.call_args
    assert kwargs.get("feature_elimination") is True
    assert kwargs.get("threshold") == 0.3


# ── compute_gene_pairs() tests ─────────────────────────────────────────────────


def test_compute_gene_pairs_marks_step_complete(adata_spatial, monkeypatch):
    ha = _setup_ha_no_network(adata_spatial, monkeypatch)
    with patch("scvi.external.harreman._analysis._compute_gene_pairs") as mock_gp:
        mock_gp.return_value = None
        ha.compute_gene_pairs()
    assert STEP_GENE_PAIRS in ha._completed_steps


def test_compute_gene_pairs_passes_layer_key(adata_spatial, monkeypatch):
    import scvi.external.harreman._analysis as _mod

    monkeypatch.setattr(_mod, "_extract_interaction_db", _make_mock_extract_db(adata_spatial))
    adata_spatial.layers["my_layer"] = adata_spatial.X.copy()
    ha = HarremanAnalysis(adata_spatial, layer_key="my_layer")
    ha.setup(compute_neighbors_on_key="spatial", n_neighbors=5)
    with patch("scvi.external.harreman._analysis._compute_gene_pairs") as mock_gp:
        mock_gp.return_value = None
        ha.compute_gene_pairs()
    _, kwargs = mock_gp.call_args
    assert kwargs.get("layer_key") == "my_layer"


# ── compute_cell_communication() tests ────────────────────────────────────────


def _ha_after_gene_pairs(adata, monkeypatch):
    ha = _setup_ha_no_network(adata, monkeypatch, cell_type_key="cell_type")
    ha._completed_steps.add(STEP_GENE_PAIRS)
    return ha


def test_ccc_standard_marks_step(adata_spatial, monkeypatch):
    ha = _ha_after_gene_pairs(adata_spatial, monkeypatch)
    with patch("scvi.external.harreman._analysis._compute_cell_communication") as mock_ccc:
        mock_ccc.return_value = None
        ha.compute_cell_communication(mode="standard")
    assert STEP_CCC in ha._completed_steps


def test_ccc_cell_type_mode_calls_ct_function(adata_spatial, monkeypatch):
    ha = _ha_after_gene_pairs(adata_spatial, monkeypatch)
    with patch("scvi.external.harreman._analysis._compute_ct_cell_communication") as mock_ct:
        mock_ct.return_value = None
        ha.compute_cell_communication(mode="cell_type")
    mock_ct.assert_called_once()


def test_ccc_standard_mode_calls_standard_function(adata_spatial, monkeypatch):
    ha = _ha_after_gene_pairs(adata_spatial, monkeypatch)
    with patch("scvi.external.harreman._analysis._compute_cell_communication") as mock_std:
        mock_std.return_value = None
        ha.compute_cell_communication(mode="standard")
    mock_std.assert_called_once()


def test_ccc_invalid_mode_raises(adata_spatial, monkeypatch):
    ha = _ha_after_gene_pairs(adata_spatial, monkeypatch)
    with pytest.raises(ValueError, match="mode"):
        ha.compute_cell_communication(mode="invalid")


# ── compute_interacting_cell_scores() + select_significant_interactions() ──────


def _ha_after_ccc(adata, monkeypatch):
    ha = _ha_after_gene_pairs(adata, monkeypatch)
    ha._completed_steps.add(STEP_CCC)
    return ha


def test_ics_standard_marks_step(adata_spatial, monkeypatch):
    ha = _ha_after_ccc(adata_spatial, monkeypatch)
    with patch("scvi.external.harreman._analysis._compute_interacting_cell_scores") as mock_ics:
        mock_ics.return_value = None
        ha.compute_interacting_cell_scores(mode="standard")
    assert STEP_ICS in ha._completed_steps


def test_ics_cell_type_mode_calls_ct_function(adata_spatial, monkeypatch):
    ha = _ha_after_ccc(adata_spatial, monkeypatch)
    with patch("scvi.external.harreman._analysis._compute_ct_interacting_cell_scores") as mock_ct:
        mock_ct.return_value = None
        ha.compute_interacting_cell_scores(mode="cell_type")
    mock_ct.assert_called_once()


def test_select_significant_marks_step(adata_spatial, monkeypatch):
    ha = _ha_after_ccc(adata_spatial, monkeypatch)
    with patch("scvi.external.harreman._analysis._select_significant_interactions") as mock_sig:
        mock_sig.return_value = None
        ha.select_significant_interactions()
    assert STEP_SIG in ha._completed_steps


def test_select_significant_passes_threshold(adata_spatial, monkeypatch):
    ha = _ha_after_ccc(adata_spatial, monkeypatch)
    with patch("scvi.external.harreman._analysis._select_significant_interactions") as mock_sig:
        mock_sig.return_value = None
        ha.select_significant_interactions(fdr_threshold=0.01)
    _, kwargs = mock_sig.call_args
    assert kwargs.get("threshold") == 0.01


def test_select_significant_defaults_to_standard_mode(adata_spatial, monkeypatch):
    ha = _ha_after_gene_pairs(adata_spatial, monkeypatch)
    with patch("scvi.external.harreman._analysis._compute_cell_communication") as mock_ccc:
        mock_ccc.return_value = None
        ha.compute_cell_communication(mode="standard")
    with patch("scvi.external.harreman._analysis._select_significant_interactions") as mock_sig:
        mock_sig.return_value = None
        ha.select_significant_interactions()
    _, kwargs = mock_sig.call_args
    assert kwargs.get("ct_aware") is False


def test_select_significant_defaults_to_cell_type_mode(adata_spatial, monkeypatch):
    ha = _ha_after_gene_pairs(adata_spatial, monkeypatch)
    with patch("scvi.external.harreman._analysis._compute_ct_cell_communication") as mock_ccc:
        mock_ccc.return_value = None
        ha.compute_cell_communication(mode="cell_type")
    with patch("scvi.external.harreman._analysis._select_significant_interactions") as mock_sig:
        mock_sig.return_value = None
        ha.select_significant_interactions()
    _, kwargs = mock_sig.call_args
    assert kwargs.get("ct_aware") is True


def test_ics_mode_mismatch_warns(adata_spatial, monkeypatch):
    ha = _ha_after_gene_pairs(adata_spatial, monkeypatch)
    with patch("scvi.external.harreman._analysis._compute_cell_communication") as mock_ccc:
        mock_ccc.return_value = None
        ha.compute_cell_communication(mode="standard")
    with patch(
        "scvi.external.harreman._analysis._compute_ct_interacting_cell_scores"
    ) as mock_ct_ics:
        mock_ct_ics.return_value = None
        with pytest.warns(UserWarning, match="does not match"):
            ha.compute_interacting_cell_scores(mode="cell_type")


def test_repr_contains_key_fields(adata_spatial):
    ha = HarremanAnalysis(adata_spatial)
    r = repr(ha)
    assert "is_set_up=False" in r
    assert "is_deconvolved=False" in r
    assert "completed_steps=" in r


def test_repr_updates_after_setup(adata_spatial, monkeypatch):
    ha = _setup_ha_no_network(adata_spatial, monkeypatch)
    r = repr(ha)
    assert "is_set_up=True" in r
    assert STEP_SETUP in r


# ── Integration tests (require model training) ────────────────────────────────
# Run with: pytest -m optional tests/external/harreman/test_harreman_analysis.py


@pytest.fixture(scope="session")
def destvi_model_and_adata():
    """Train minimal DestVI model on synthetic data."""
    from scvi.data import synthetic_iid
    from scvi.model import CondSCVI, DestVI

    n_labels = 3
    dataset = synthetic_iid(n_labels=n_labels, n_genes=50)
    dataset.obs["overclustering_vamp"] = list(range(dataset.n_obs))
    CondSCVI.setup_anndata(dataset, labels_key="labels")
    sc_model = CondSCVI(dataset, n_latent=2, n_layers=1, prior="mog", num_classes_mog=10)
    sc_model.train(1, train_size=1)
    del dataset.obs["overclustering_vamp"]
    DestVI.setup_anndata(dataset, layer=None)
    spatial_model = DestVI.from_rna_model(dataset, sc_model, amortization="latent")
    spatial_model.train(max_epochs=1)
    dataset.obsm["spatial"] = np.random.default_rng(0).random((dataset.n_obs, 2)) * 100
    return spatial_model, dataset


@pytest.fixture(scope="session")
def resolvi_model_and_adata():
    """Train minimal RESOLVI model on synthetic spatial data."""
    from scvi.data import synthetic_iid
    from scvi.external import RESOLVI

    adata = synthetic_iid(generate_coordinates=True, n_regions=5, n_genes=50)
    adata.obsm["X_spatial"] = adata.obsm["coordinates"]
    RESOLVI.setup_anndata(adata)
    model = RESOLVI(adata)
    model.train(max_epochs=2)
    return model, adata


@pytest.fixture(scope="session")
def scviva_model_and_adata():
    """Train minimal SCVIVA model on synthetic spatial data."""
    from scvi.data import synthetic_iid
    from scvi.external import SCVIVA

    adata = synthetic_iid(
        batch_size=64,
        n_genes=50,
        n_proteins=0,
        n_regions=0,
        n_batches=2,
        n_labels=3,
        generate_coordinates=True,
        sparse_format=None,
    )
    adata.layers["counts"] = adata.X.copy()
    adata.obsm["qz1_m"] = np.random.default_rng(1).normal(size=(adata.n_obs, 10))
    SCVIVA.preprocessing_anndata(
        adata,
        k_nn=5,
        sample_key="batch",
        labels_key="labels",
        cell_coordinates_key="coordinates",
        expression_embedding_key="qz1_m",
        expression_embedding_niche_key="qz1_m_niche_ct",
        niche_composition_key="neighborhood_composition",
        niche_indexes_key="niche_indexes",
        niche_distances_key="niche_distances",
    )
    SCVIVA.setup_anndata(
        adata,
        layer="counts",
        batch_key="batch",
        sample_key="batch",
        labels_key="labels",
        cell_coordinates_key="coordinates",
        expression_embedding_key="qz1_m",
        expression_embedding_niche_key="qz1_m_niche_ct",
        niche_composition_key="neighborhood_composition",
        niche_indexes_key="niche_indexes",
        niche_distances_key="niche_distances",
    )
    model = SCVIVA(adata, prior_mixture=False, semisupervised=False)
    model.train(max_epochs=2, accelerator="cpu")
    return model, adata


@pytest.mark.optional
def test_integration_destvi_attaches_ct_layers(destvi_model_and_adata):
    """DestVI: verify each cell-type layer is attached to adata."""
    model, adata = destvi_model_and_adata
    ha = HarremanAnalysis(adata.copy(), model=model)
    proportions = model.get_proportions()
    for ct in proportions.columns:
        assert ct in ha.adata.layers, f"Missing layer for cell type '{ct}'"


@pytest.mark.optional
def test_integration_destvi_is_deconvolved(destvi_model_and_adata):
    model, adata = destvi_model_and_adata
    ha = HarremanAnalysis(adata.copy(), model=model)
    assert ha.is_deconvolved is True


@pytest.mark.optional
def test_integration_destvi_layer_shape(destvi_model_and_adata):
    """Cell-type layers must have same shape as adata.X."""
    model, adata = destvi_model_and_adata
    adata_copy = adata.copy()
    ha = HarremanAnalysis(adata_copy, model=model)
    proportions = model.get_proportions()
    for ct in proportions.columns:
        assert ha.adata.layers[ct].shape == adata.shape, f"Layer '{ct}' shape mismatch"


@pytest.mark.optional
def test_integration_resolvi_attaches_denoised_layer(resolvi_model_and_adata):
    """RESOLVI: verify denoised layer attached."""
    from scvi.external.harreman._constants import HARREMAN_DENOISED_LAYER

    model, adata = resolvi_model_and_adata
    ha = HarremanAnalysis(adata.copy(), model=model)
    assert HARREMAN_DENOISED_LAYER in ha.adata.layers


@pytest.mark.optional
def test_integration_resolvi_not_deconvolved(resolvi_model_and_adata):
    model, adata = resolvi_model_and_adata
    ha = HarremanAnalysis(adata.copy(), model=model)
    assert ha.is_deconvolved is False


@pytest.mark.optional
def test_integration_resolvi_denoised_shape(resolvi_model_and_adata):
    from scvi.external.harreman._constants import HARREMAN_DENOISED_LAYER

    model, adata = resolvi_model_and_adata
    adata_copy = adata.copy()
    ha = HarremanAnalysis(adata_copy, model=model)
    assert ha.adata.layers[HARREMAN_DENOISED_LAYER].shape == adata.shape


@pytest.mark.optional
def test_integration_scviva_attaches_latent(scviva_model_and_adata):
    """SCVIVA: verify latent representation attached to obsm."""
    from scvi.external.harreman._constants import HARREMAN_LATENT_OBSM

    model, adata = scviva_model_and_adata
    ha = HarremanAnalysis(adata.copy(), model=model)
    assert HARREMAN_LATENT_OBSM in ha.adata.obsm


@pytest.mark.optional
def test_integration_scviva_not_deconvolved(scviva_model_and_adata):
    model, adata = scviva_model_and_adata
    ha = HarremanAnalysis(adata.copy(), model=model)
    assert ha.is_deconvolved is False


@pytest.mark.optional
def test_integration_scviva_latent_shape(scviva_model_and_adata):
    from scvi.external.harreman._constants import HARREMAN_LATENT_OBSM

    model, adata = scviva_model_and_adata
    adata_copy = adata.copy()
    ha = HarremanAnalysis(adata_copy, model=model)
    latent = ha.adata.obsm[HARREMAN_LATENT_OBSM]
    assert latent.shape[0] == adata.n_obs
    assert latent.ndim == 2


# ── Accessor tests ─────────────────────────────────────────────────────────────


def test_accessors_exist(adata_spatial):
    ha = HarremanAnalysis(adata_spatial)
    assert hasattr(ha, "hs")
    assert hasattr(ha, "tl")
    assert hasattr(ha, "pl")


def test_hs_accessor_methods(adata_spatial):
    ha = HarremanAnalysis(adata_spatial)
    for method in [
        "compute_local_autocorrelation",
        "compute_local_correlation",
        "create_modules",
        "calculate_module_scores",
        "calculate_super_module_scores",
        "compute_top_scoring_modules",
    ]:
        assert callable(getattr(ha.hs, method)), f"ha.hs.{method} not callable"


def test_tl_accessor_methods(adata_spatial):
    ha = HarremanAnalysis(adata_spatial)
    for method in [
        "compute_knn_graph",
        "compute_interaction_module_correlation",
        "select_significant_interactions",
    ]:
        assert callable(getattr(ha.tl, method)), f"ha.tl.{method} not callable"


def test_pl_accessor_methods(adata_spatial):
    ha = HarremanAnalysis(adata_spatial)
    for method in [
        "local_correlation_plot",
        "average_local_correlation_plot",
        "module_score_correlation_plot",
        "plot_interacting_cell_scores",
        "plot_ct_interacting_cell_scores",
        "plot_interaction_module_correlation",
    ]:
        assert callable(getattr(ha.pl, method)), f"ha.pl.{method} not callable"


def test_accessor_resolve_defaults_to_ha_adata(adata_spatial):
    ha = HarremanAnalysis(adata_spatial)
    assert ha.hs._resolve(None) is adata_spatial
    assert ha.tl._resolve(None) is adata_spatial
    assert ha.pl._resolve(None) is adata_spatial


def test_accessor_resolve_accepts_override(adata_spatial):
    import anndata as ad
    other = ad.AnnData(adata_spatial.X.copy())
    ha = HarremanAnalysis(adata_spatial)
    assert ha.hs._resolve(other) is other
    assert ha.tl._resolve(other) is other
    assert ha.pl._resolve(other) is other


def test_tl_compute_knn_graph_via_accessor(adata_spatial):
    """ha.tl.compute_knn_graph runs and writes to obsp['distances']."""
    ha = HarremanAnalysis(adata_spatial)
    ha.tl.compute_knn_graph(compute_neighbors_on_key="spatial", n_neighbors=5)
    assert "distances" in ha.adata.obsp


# ── GPU / rapids_singlecell integration tests ─────────────────────────────────


def test_compute_neighbors_use_gpu_false_skips_gpu(adata_spatial):
    """use_gpu=False must not call _gpu_neighbors."""
    from scvi.external.harreman.tools import knn as knn_mod

    with patch.object(knn_mod, "_gpu_neighbors", wraps=knn_mod._gpu_neighbors) as mock_gpu:
        knn_mod.compute_neighbors(
            adata_spatial,
            compute_neighbors_on_key="spatial",
            n_neighbors=5,
            use_gpu=False,
        )
    mock_gpu.assert_not_called()
    assert "distances" in adata_spatial.obsp


def test_compute_neighbors_use_gpu_none_falls_back_when_gpu_unavailable(adata_spatial):
    """use_gpu=None falls back to sklearn when _gpu_neighbors returns None."""
    from scvi.external.harreman.tools import knn as knn_mod

    with patch.object(knn_mod, "_gpu_neighbors", return_value=None) as mock_gpu:
        knn_mod.compute_neighbors(
            adata_spatial,
            compute_neighbors_on_key="spatial",
            n_neighbors=5,
            use_gpu=None,
        )
    mock_gpu.assert_called_once()
    assert "distances" in adata_spatial.obsp


def test_compute_neighbors_use_gpu_true_raises_when_gpu_unavailable(adata_spatial):
    """use_gpu=True raises RuntimeError when _gpu_neighbors returns None."""
    from scvi.external.harreman.tools import knn as knn_mod

    with patch.object(knn_mod, "_gpu_neighbors", return_value=None):
        with pytest.raises(RuntimeError, match="GPU neighbor computation failed"):
            knn_mod.compute_neighbors(
                adata_spatial,
                compute_neighbors_on_key="spatial",
                n_neighbors=5,
                use_gpu=True,
            )


def test_compute_neighbors_uses_gpu_result_when_available(adata_spatial):
    """When _gpu_neighbors returns a matrix, it must be used without sklearn fallback."""
    import numpy as np
    from scipy.sparse import csr_matrix

    from scvi.external.harreman.tools import knn as knn_mod

    n = adata_spatial.n_obs
    # Build a minimal sparse distance matrix: each cell connected to next
    row = np.arange(n - 1)
    col = np.arange(1, n)
    data = np.ones(n - 1, dtype=np.float32)
    fake_dist = csr_matrix((data, (row, col)), shape=(n, n))

    with patch.object(knn_mod, "_gpu_neighbors", return_value=fake_dist) as mock_gpu:
        with patch.object(knn_mod, "NearestNeighbors") as mock_nn:
            knn_mod.compute_neighbors(
                adata_spatial,
                compute_neighbors_on_key="spatial",
                n_neighbors=5,
                use_gpu=None,
            )
    mock_gpu.assert_called_once()
    mock_nn.assert_not_called()  # sklearn path must be bypassed
    assert "distances" in adata_spatial.obsp


def test_setup_use_gpu_param_forwarded(adata_spatial, monkeypatch):
    """use_gpu kwarg is forwarded from ha.setup() to compute_knn_graph."""
    import scvi.external.harreman._analysis as _mod

    monkeypatch.setattr(_mod, "_extract_interaction_db", _make_mock_extract_db(adata_spatial))
    with patch.object(_mod, "_compute_knn_graph") as mock_knn:
        mock_knn.return_value = None
        ha = HarremanAnalysis(adata_spatial)
        ha.setup(compute_neighbors_on_key="spatial", n_neighbors=5, use_gpu=False)
    _, kwargs = mock_knn.call_args
    assert kwargs.get("use_gpu") is False
