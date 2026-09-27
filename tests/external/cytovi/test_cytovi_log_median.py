import warnings

import anndata as ad
import numpy as np
import pandas as pd
import pytest

from scvi.external import CYTOVI
from scvi.external.cytovi._utils import log_median


@pytest.mark.parametrize("shape", [(3, 5), (4, 6)])
@pytest.mark.parametrize("axis", [0, 1, -1, None, (0, 1)])
def test_log_median_matches_probability_space_statistic(shape, axis):
    values = np.random.default_rng(42).uniform(-20, 20, size=shape)
    original = values.copy()
    expected = np.log(np.median(np.exp(values), axis=axis))
    np.testing.assert_allclose(log_median(values, axis=axis), expected, atol=1e-12)
    np.testing.assert_array_equal(values, original)


@pytest.mark.parametrize("dtype", [np.float32, np.float64])
@pytest.mark.parametrize("count", [1, 5, 6])
@pytest.mark.parametrize("offset", [-2000.0, 2000.0])
def test_log_median_is_stable_under_large_shifts(dtype, count, offset):
    values = np.arange(2 * count, dtype=dtype).reshape(2, count)
    expected = np.log(np.median(np.exp(values.astype(float)), axis=1)) + offset
    with np.errstate(over="raise", divide="raise", invalid="raise"):
        actual = log_median(values + offset)
    np.testing.assert_allclose(actual, expected, rtol=1e-6)


def test_even_log_median_keeps_probability_space_average():
    values = np.log([[1.0, 9.0]])
    np.testing.assert_allclose(log_median(values), np.log([5.0]))
    assert not np.allclose(log_median(values), np.median(values, axis=1))


def test_log_median_with_widely_separated_values():
    values = [[-2000.0, -2000.0, 2000.0], [-2000.0, 2000.0, 2000.0]]
    np.testing.assert_allclose(log_median(values), [-2000.0, 2000.0])
    np.testing.assert_allclose(log_median([[-2000.0, 2000.0]]), [2000.0 - np.log(2)])


@pytest.mark.parametrize(
    ("values", "expected"),
    [
        ([-np.inf, -np.inf], -np.inf),
        ([np.inf, np.inf], np.inf),
        ([-np.inf, np.inf], np.inf),
        ([-np.inf, 0.0], -np.log(2)),
        ([-np.inf, 0.0, np.inf], 0.0),
        ([0.0, 1.0, np.nan], np.nan),
        ([0.0, 1.0, 2.0, np.nan], np.nan),
    ],
)
def test_log_median_preserves_nonfinite_input_semantics(values, expected):
    np.testing.assert_allclose(log_median([values]), [expected], equal_nan=True)


@pytest.mark.parametrize(("shape", "axis"), [((0, 3), 0), ((0, 3), 1), ((3, 0), 1)])
def test_log_median_preserves_empty_reduction(shape, axis):
    values = np.empty(shape)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", RuntimeWarning)
        expected = np.log(np.median(np.exp(values), axis=axis))
        actual = log_median(values, axis=axis)
    np.testing.assert_allclose(actual, expected, equal_nan=True)


@pytest.mark.parametrize("samples_per_group", [1, 2, 3])
def test_da_remains_finite_for_extreme_log_probabilities(monkeypatch, samples_per_group):
    samples = [f"sample_{i}" for i in range(2 * samples_per_group)]
    adata = ad.AnnData(np.ones((len(samples), 2), dtype=np.float32))
    adata.obs["sample"] = samples
    adata.obs["condition"] = ["case"] * samples_per_group + ["control"] * samples_per_group
    CYTOVI.setup_anndata(adata, sample_key="sample")
    model = CYTOVI(adata, n_latent=2)
    scores = pd.DataFrame(
        np.tile([-900.0] * samples_per_group + [-880.0] * samples_per_group, (adata.n_obs, 1)),
        index=adata.obs_names,
        columns=samples,
    )
    monkeypatch.setattr(model, "get_sample_logprobs", lambda *args, **kwargs: scores)
    with np.errstate(over="raise", divide="raise", invalid="raise"):
        actual = model.differential_abundance(groupby="condition")
    np.testing.assert_allclose(actual["DA_case"], -20.0)
    np.testing.assert_allclose(actual["DA_control"], 20.0)
