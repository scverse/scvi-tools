from __future__ import annotations

from typing import TYPE_CHECKING

import numpy as np
import pandas as pd
import pytest
import torch.distributions as dist
from mudata import MuData

from scvi.data import synthetic_iid
from scvi.external import TOTALANVI
from scvi.model import MULTIVI, SCANVI, SCVI, TOTALVI, AmortizedLDA

if TYPE_CHECKING:
    from anndata import AnnData

    from scvi.model.base import VAEMixin


@pytest.fixture(scope="session")
def adata():
    adata = synthetic_iid(batch_size=500)
    adata.obs["sample"] = np.random.choice(10, size=adata.n_obs)
    adata.obs["sample_str"] = [chr(i + ord("a")) for i in adata.obs["sample"]]
    return adata


@pytest.fixture(scope="session")
def mdata():
    adata = synthetic_iid(batch_size=500)
    protein_adata = synthetic_iid(batch_size=500)
    protein_adata.var_names = protein_adata.uns["protein_names"]
    mdata = MuData({"rna": adata, "protein": protein_adata})
    mdata.obs["sample"] = np.random.choice(10, size=mdata.n_obs)
    mdata.obs["sample_str"] = [chr(i + ord("a")) for i in mdata.obs["sample"]]
    return mdata


@pytest.fixture(
    scope="session",
    params=[SCVI, SCANVI, TOTALVI, MULTIVI, TOTALANVI],
)
def model(request, adata, mdata):
    model_cls = request.param
    model_kwargs = {}

    if model_cls is SCVI:
        model_cls.setup_anndata(adata=adata, batch_key="batch")
    elif model_cls is SCANVI:
        model_cls.setup_anndata(
            adata=adata, labels_key="labels", unlabeled_category="NA", batch_key="batch"
        )
    elif model_cls in (TOTALVI, TOTALANVI, MULTIVI):
        adata = mdata
        setup_kwargs = {}
        modalities = {"rna_layer": "rna", "batch_key": "rna", "protein_layer": "protein"}
        if model_cls is TOTALANVI:
            setup_kwargs = {"labels_key": "labels", "unlabeled_category": "label_0"}
            modalities["labels_key"] = "rna"
        if model_cls in (TOTALVI, TOTALANVI):
            model_kwargs["empirical_protein_background_prior"] = False
        model_cls.setup_mudata(
            mdata=adata,
            batch_key="batch",
            modalities=modalities,
            **setup_kwargs,
        )
    else:
        model_cls.setup_anndata(adata=adata)

    model_inst = model_cls(adata, **model_kwargs)
    model_inst.train(max_epochs=1, train_size=0.5)
    return model_inst


@pytest.mark.parametrize(
    "ap_kwargs",
    [
        {},
        {"indices": None},
        {"indices": np.array([])},
        {"indices": np.random.choice(500, 100, replace=False)},
        {"indices": np.arange(150), "dof": None},
        {"indices": np.arange(150), "dof": 5},
    ],
)
def test_get_aggregated_posterior(model: VAEMixin, adata: AnnData, mdata: MuData, ap_kwargs):
    if isinstance(model.adata, MuData):
        adata = mdata

    ap = model.get_aggregated_posterior(adata, **ap_kwargs)
    assert isinstance(ap, dist.Distribution)

    subset_indices = np.random.choice(adata.n_obs, adata.n_obs // 2, replace=False)
    adata_subset = adata[subset_indices, :].copy()
    ap = model.get_aggregated_posterior(adata_subset, **ap_kwargs)
    assert isinstance(ap, dist.Distribution)


@pytest.mark.parametrize(
    "da_kwargs",
    [
        {"sample_key": "sample_str", "num_cells_posterior": 100, "dof": 3},
        {"sample_key": "sample_str", "num_cells_posterior": 5000, "dof": 3},
        {"sample_key": "sample_str", "num_cells_posterior": None, "dof": 3},
        {"sample_key": "sample_str", "dof": 3},
        {"sample_key": "sample", "num_cells_posterior": 100, "dof": 3},
        {"sample_key": None, "num_cells_posterior": 100, "dof": 3},
        {"sample_key": "sample_str", "num_cells_posterior": 100},
        {"sample_key": "sample_str", "num_cells_posterior": 100, "dof": None},
    ],
)
def test_differential_abundance(model: VAEMixin, adata: AnnData, mdata: MuData, da_kwargs):
    if isinstance(model.adata, MuData):
        adata = mdata
    if da_kwargs["sample_key"] is None:
        with pytest.raises(KeyError):
            model.differential_abundance(adata, **da_kwargs)
    else:
        model.differential_abundance(adata, **da_kwargs)
        assert isinstance(adata.obsm["da_log_probs"], pd.DataFrame)

        subset_indices = np.random.choice(adata.n_obs, adata.n_obs // 2, replace=False)
        adata_subset = adata[subset_indices, :].copy()
        model.differential_abundance(adata, adata_subset, **da_kwargs)
        assert isinstance(adata_subset.obsm["da_log_probs"], pd.DataFrame)


def test_differential_abundance_not_implemented(adata):
    AmortizedLDA.setup_anndata(adata)
    model = AmortizedLDA(adata, n_topics=5)
    model.train(max_epochs=1)
    with pytest.raises(
        NotImplementedError, match="differential_abundance is not implemented for AmortizedLDA"
    ):
        model.differential_abundance(adata)
