from __future__ import annotations

from typing import TYPE_CHECKING

import numpy as np
import pandas as pd
import pytest
import torch.distributions as dist

from scvi.data import synthetic_iid
from scvi.model import SCVI
from scvi.model.base._da_pytorch import differential_abundance, get_aggregated_posterior

if TYPE_CHECKING:
    from anndata import AnnData


@pytest.fixture(scope="session")
def adata():
    adata = synthetic_iid(batch_size=500)
    adata.obs["sample"] = np.random.choice(10, size=adata.n_obs)
    adata.obs["sample_str"] = [chr(i + ord("a")) for i in adata.obs["sample"]]
    return adata


@pytest.fixture(scope="session")
def model(adata):
    SCVI.setup_anndata(adata=adata, batch_key="batch")
    model = SCVI(adata)
    model.train(max_epochs=1, train_size=0.5)
    return model


@pytest.mark.parametrize(
    "ap_kwargs",
    [
        {},
        {"indices": None},
        {"indices": []},
        {"indices": np.random.choice(1500, 100, replace=False)},
        {"indices": list(range(150)), "dof": None},
        {"indices": list(range(150)), "dof": 5},
    ],
)
def test_get_aggregated_posterior(model: SCVI, adata: AnnData, ap_kwargs):
    ap = get_aggregated_posterior(model, adata, **ap_kwargs)
    assert isinstance(ap, dist.Distribution)

    subset_indices = np.random.choice(adata.n_obs, adata.n_obs // 2, replace=False)
    adata_subset = adata[subset_indices, :].copy()
    ap = get_aggregated_posterior(model, adata_subset, **ap_kwargs)
    assert isinstance(ap, dist.Distribution)


@pytest.mark.parametrize(
    "da_kwargs",
    [
        {"adata": adata, "sample_key": "sample_str", "num_cells_posterior": 1000, "dof": 3},
    ],
)
def test_differential_abundance(model: SCVI, adata: AnnData, da_kwargs):
    differential_abundance(model, adata, **da_kwargs)
    assert isinstance(adata.obsm["da_log_probs"], pd.DataFrame)

    subset_indices = np.random.choice(adata.n_obs, adata.n_obs // 2, replace=False)
    adata_subset = adata[subset_indices, :].copy()
    differential_abundance(model, adata_subset, **da_kwargs)
    assert isinstance(adata_subset.obsm["da_log_probs"], pd.DataFrame)
