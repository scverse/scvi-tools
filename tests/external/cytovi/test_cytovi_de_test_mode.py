import anndata as ad
import numpy as np
import pytest

from scvi.external import CYTOVI


@pytest.fixture
def model():
    adata = ad.AnnData(np.random.default_rng(42).uniform(0.1, 0.9, (12, 4)).astype(np.float32))
    adata.obs_names = [f"cell_{i}" for i in range(12)]
    adata.obs["sample"] = ["a"] * 6 + ["b"] * 6
    CYTOVI.setup_anndata(adata, sample_key="sample")
    model = CYTOVI(adata, n_latent=2)
    model.is_trained_ = True
    model.module.eval()
    return model


@pytest.mark.parametrize("lfc_clipping", [False, True])
@pytest.mark.parametrize("test_mode", ["two", "three"])
def test_de_forwards_test_mode_independently_of_clipping(
    model, monkeypatch, lfc_clipping, test_mode
):
    from scvi.external.cytovi import _model

    received = {}
    sentinel = object()

    def capture_de(*args, **kwargs):
        received.update(kwargs)
        return sentinel

    monkeypatch.setattr(_model, "_de_core", capture_de)
    actual = model.differential_expression(
        groupby="sample", test_mode=test_mode, lfc_clipping=lfc_clipping, balance_samples=False
    )
    assert actual is sentinel
    assert received["test_mode"] == test_mode
    assert ("change_fn" in received) == lfc_clipping


def test_de_default_test_mode_is_two_without_clipping(model, monkeypatch):
    from scvi.external.cytovi import _model

    received = {}
    monkeypatch.setattr(_model, "_de_core", lambda *args, **kwargs: received.update(kwargs))
    model.differential_expression(groupby="sample", lfc_clipping=False, balance_samples=False)
    assert received["test_mode"] == "two"
