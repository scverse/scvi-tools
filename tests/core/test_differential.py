from functools import partial

import numpy as np
import pytest

from scvi.data import synthetic_iid
from scvi.model import SCVI
from scvi.model.base._differential import (
    DifferentialComputation,
    estimate_delta,
    estimate_pseudocounts_offset,
)
from scvi.model.base._utils import _prepare_obs


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
    expected_range = (delta >= 0.4 * 3) and (delta <= 6)
    if not expected_range:
        raise ValueError("The effect-size threshold was not properly estimated.")

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
    expected_off_range = offset <= 1e-6
    if not expected_off_range:
        raise ValueError("The pseudocount offset was not properly estimated.")


def test_differential_computation(save_path):

    n_latent = 5
    adata = synthetic_iid()
    SCVI.setup_anndata(
        adata,
        batch_key="batch",
        labels_key="labels",
    )
    model = SCVI(adata, n_latent=n_latent)
    model.train(1)

    model_fn = partial(model.get_normalized_expression, return_numpy=True)
    dc = DifferentialComputation(model_fn, model.adata_manager)

    cell_idx1 = np.asarray(adata.obs.labels == "label_1")
    cell_idx2 = ~cell_idx1

    dc.get_bayes_factors(cell_idx1, cell_idx2, mode="vanilla", use_permutation=True)
    res = dc.get_bayes_factors(
        cell_idx1, cell_idx2, mode="change", use_permutation=False
    )
    assert (res["delta"] == 0.5) and (res["pseudocounts"] == 0.0)
    res = dc.get_bayes_factors(
        cell_idx1, cell_idx2, mode="change", use_permutation=False, delta=None
    )
    dc.get_bayes_factors(
        cell_idx1,
        cell_idx2,
        mode="change",
        use_permutation=False,
        delta=None,
        pseudocounts=None,
    )
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

    # Test query features
    obs_col, group1, _, = _prepare_obs(
        idx1="(labels == 'label_1') & (batch == 'batch_1')", idx2=None, adata=adata
    )
    assert (obs_col == group1).sum() == adata.obs.loc[
        lambda x: (x.labels == "label_1") & (x.batch == "batch_1")
    ].shape[0]
    model.differential_expression(
        idx1="labels == 'label_1'",
    )
    model.differential_expression(
        idx1="labels == 'label_1'", idx2="(labels == 'label_2') & (batch == 'batch_1')"
    )

    # test that ints as group work
    a = synthetic_iid()
    SCVI.setup_anndata(
        a,
        batch_key="batch",
        labels_key="labels",
    )
    a.obs["test"] = [0] * 200 + [1] * 200
    model = SCVI(a)
    model.differential_expression(groupby="test", group1=0)

    # test that string but not as categorical work
    a = synthetic_iid()
    SCVI.setup_anndata(
        a,
        batch_key="batch",
        labels_key="labels",
    )
    a.obs["test"] = ["0"] * 200 + ["1"] * 200
    model = SCVI(a)
    model.differential_expression(groupby="test", group1="0")
