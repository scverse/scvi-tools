import anndata as ad
import numpy as np
import pandas as pd
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


@pytest.fixture
def log_probs(model, monkeypatch):
    # Deliberately reverse sample columns relative to their first appearance in obs.
    values = pd.DataFrame(
        np.column_stack([np.arange(12) - 20.0, np.arange(12) - 10.0]),
        index=model.adata.obs_names,
        columns=["b", "a"],
    )
    monkeypatch.setattr(model, "get_sample_logprobs", lambda *args, **kwargs: values)
    model.adata.obs["condition"] = ["case"] * 6 + ["control"] * 6
    return values


def test_da_preserves_cell_names_and_aligns_sample_conditions(model, log_probs):
    actual = model.differential_abundance(groupby="condition")
    assert actual.index.equals(log_probs.index)
    np.testing.assert_allclose(actual["DA_case"], 10.0)
    np.testing.assert_allclose(actual["DA_control"], -10.0)


def test_da_can_group_by_sample(model, log_probs):
    actual = model.differential_abundance(groupby="sample")
    assert actual.index.equals(log_probs.index)
    np.testing.assert_allclose(actual["DA_a"], 10.0)
    np.testing.assert_allclose(actual["DA_b"], -10.0)


@pytest.mark.parametrize("invalid", ["mixed", "missing", "missing_sample", "single"])
def test_da_rejects_invalid_sample_conditions(model, log_probs, invalid):
    if invalid == "mixed":
        model.adata.obs.loc["cell_0", "condition"] = "control"
        match = "one condition per sample"
    elif invalid == "missing":
        model.adata.obs.loc["cell_0", "condition"] = None
        match = "missing"
    elif invalid == "missing_sample":
        model.adata.obs.loc["cell_0", "sample"] = None
        match = "missing"
    else:
        model.adata.obs["condition"] = "case"
        match = "at least two"
    with pytest.raises(ValueError, match=match):
        model.differential_abundance(groupby="condition")


def test_da_log_probabilities_can_be_returned_without_conditions(model, log_probs):
    actual = model.differential_abundance(groupby="unused", return_log_probs=True)
    pd.testing.assert_frame_equal(actual, log_probs)
