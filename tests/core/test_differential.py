from functools import partial

import numpy as np
import pytest

from scvi.data import synthetic_iid
from scvi.model import SCVI
from scvi.utils import DifferentialComputation
from scvi.utils._differential import estimate_delta, estimate_pseudocounts_offset


def test_features():
    a = np.random.randn(
        100,
    )
    b = 3 + np.random.randn(
        100,
    )
    c = -3 + np.random.randn(
        100,
    )
    alls = np.concatenate([a, b, c])
    delta = estimate_delta(alls)
    assert delta >= 0.4 * 3
    assert delta <= 6

    scales_a = np.random.rand(100, 50)
    where_zero_a = np.zeros(50, dtype=bool)
    where_zero_a[:10] = True
    scales_a[:, :10] = 1e-6

    scales_b = np.random.rand(100, 50)
    where_zero_b = np.zeros(50, dtype=bool)
    where_zero_b[-10:] = True
    scales_b[:, -10:] = 1e-7
    offset = estimate_pseudocounts_offset(
        scales_a=scales_a,
        scales_b=scales_b,
        where_zero_a=where_zero_a,
        where_zero_b=where_zero_b,
    )
    assert offset <= 1e-6


def test_differential_computation(save_path):

    n_latent = 5
    adata = synthetic_iid()
    model = SCVI(adata, n_latent=n_latent)
    model.train(1)

    model_fn = partial(model.get_normalized_expression, return_numpy=True)
    dc = DifferentialComputation(model_fn, adata)

    cell_idx1 = np.asarray(adata.obs.labels == "label_1")
    cell_idx2 = ~cell_idx1

    dc.get_bayes_factors(cell_idx1, cell_idx2, mode="vanilla", use_permutation=True)
    dc.get_bayes_factors(cell_idx1, cell_idx2, mode="change", use_permutation=False)
    dc.get_bayes_factors(cell_idx1, cell_idx2, mode="change", cred_interval_lvls=[0.75])

    delta = 0.5

    def change_fn_test(x, y):
        return x - y

    def m1_domain_fn_test(samples):
        return np.abs(samples) >= delta

    dc.get_bayes_factors(
        cell_idx1,
        cell_idx2,
        mode="change",
        m1_domain_fn=m1_domain_fn_test,
        change_fn=change_fn_test,
    )

    # should fail if just one batch
    with pytest.raises(ValueError):
        model.differential_expression(adata[:20], groupby="batch")

    # test view
    model.differential_expression(
        adata[adata.obs["labels"] == "label_1"], groupby="batch"
    )

    # test that ints as group work
    a = synthetic_iid()
    a.obs["test"] = [0] * 200 + [1] * 200
    model = SCVI(a)
    model.differential_expression(groupby="test", group1=0)

    # test that string but not as categorical work
    a = synthetic_iid()
    a.obs["test"] = ["0"] * 200 + ["1"] * 200
    model = SCVI(a)
    model.differential_expression(groupby="test", group1="0")
