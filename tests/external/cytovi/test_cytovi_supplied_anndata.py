import anndata as ad
import numpy as np
import pytest
import torch

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


@pytest.mark.parametrize("dof", [None, 3.0])
def test_aggregated_posterior_defaults_to_supplied_cells(model, dof):
    adata = model.adata[[9, 1, 7, 3]].copy()
    expected = model.get_aggregated_posterior(adata, indices=np.arange(adata.n_obs), dof=dof)
    actual = model.get_aggregated_posterior(adata, dof=dof)
    torch.testing.assert_close(
        actual.component_distribution.loc, expected.component_distribution.loc
    )
    torch.testing.assert_close(
        actual.component_distribution.scale, expected.component_distribution.scale
    )


@pytest.mark.parametrize("rows", [[9, 1, 7, 3], list(reversed(range(12)))])
@pytest.mark.parametrize("dof", [None, 3.0])
def test_sample_logprobs_use_supplied_cells_and_order(model, rows, dof):
    adata = model.adata[rows].copy()
    latent = model.get_latent_representation(adata=adata, give_mean=True)
    expected = []
    for sample in adata.obs["sample"].unique():
        indices = np.flatnonzero(adata.obs["sample"] == sample)
        posterior = model.get_aggregated_posterior(adata, indices=indices, dof=dof)
        expected.append(posterior.log_prob(torch.as_tensor(latent)).sum(-1).numpy())

    actual = model.get_sample_logprobs(adata, batch_size=3, dof=dof)

    assert actual.index.tolist() == adata.obs_names.tolist()
    assert actual.columns.tolist() == adata.obs["sample"].unique().tolist()
    np.testing.assert_allclose(actual.to_numpy(), np.column_stack(expected), rtol=1e-5, atol=1e-6)
