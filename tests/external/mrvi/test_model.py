from __future__ import annotations

import os
from typing import TYPE_CHECKING

import numpy as np
import pytest

from scvi.data import synthetic_iid
from scvi.external import MRVI

if TYPE_CHECKING:
    from typing import Any

    from anndata import AnnData


@pytest.fixture(scope="session")
def adata():
    adata = synthetic_iid()
    adata.obs.index.name = "cell_id"
    adata.obs["sample"] = np.random.choice(15, size=adata.shape[0])
    adata.obs["sample_str"] = [chr(i + ord("a")) for i in adata.obs["sample"]]
    meta1 = np.random.randint(0, 2, size=15)
    adata.obs["meta1"] = meta1[adata.obs["sample"].values]
    meta2 = np.random.randn(15)
    adata.obs["meta2"] = meta2[adata.obs["sample"].values]
    adata.obs["cont_cov"] = np.random.normal(0, 1, size=adata.shape[0])
    adata.obs["meta1_cat"] = "CAT_" + adata.obs["meta1"].astype(str)
    adata.obs["meta1_cat"] = adata.obs["meta1_cat"].astype("category")
    adata.obs.loc[:, "disjoint_batch"] = (adata.obs.loc[:, "sample"] <= 6).replace(
        {True: "batch_0", False: "batch_1"}
    )
    adata.obs["dummy_batch"] = 1
    return adata


@pytest.fixture(scope="session")
def model(adata: AnnData):
    MRVI.setup_anndata(adata, sample_key="sample_str", batch_key="batch")
    model = MRVI(adata)
    model.train(max_steps=2, train_size=0.5)

    return model


@pytest.mark.optional
def test_mrvi(model: MRVI, adata: AnnData, save_path: str):
    model.get_local_sample_distances()
    model.get_local_sample_distances(normalize_distances=True)
    model.get_latent_representation(give_z=False)
    model.get_latent_representation(give_z=True)

    model_path = os.path.join(save_path, "mrvi_model")
    model.save(model_path, save_anndata=False, overwrite=True)
    model = MRVI.load(model_path, adata=adata)


@pytest.mark.optional
@pytest.mark.parametrize(
    ("setup_kwargs", "de_kwargs"),
    [
        (
            {"sample_key": "sample_str", "batch_key": "batch"},
            [
                {
                    "sample_cov_keys": ["meta1_cat", "meta2", "cont_cov"],
                    "store_lfc": True,
                    "add_batch_specific_offsets": True,
                },
                {
                    "sample_cov_keys": ["meta1_cat", "meta2", "cont_cov"],
                    "store_lfc": True,
                    "add_batch_specific_offsets": True,
                    "filter_inadmissible_samples": True,
                },
                {
                    "sample_cov_keys": ["meta1_cat", "meta2", "cont_cov"],
                    "store_lfc": True,
                    "add_batch_specific_offsets": False,
                },
            ],
        ),
        (
            {"sample_key": "sample_str", "batch_key": "dummy_batch"},
            [
                {
                    "sample_cov_keys": ["meta1_cat", "meta2", "cont_cov"],
                    "store_lfc": True,
                },
                {
                    "sample_cov_keys": ["meta1_cat", "meta2", "cont_cov"],
                    "store_lfc": True,
                    "lambd": 1e-1,
                },
                {
                    "sample_cov_keys": ["meta1_cat", "meta2", "cont_cov"],
                    "store_lfc": True,
                    "filter_inadmissible_samples": True,
                },
            ],
        ),
    ],
)
def test_mrvi_de(model: MRVI, setup_kwargs: dict[str, Any], de_kwargs: dict[str, Any]):
    for de_kwarg in de_kwargs:
        model.differential_expression(**de_kwarg)


@pytest.mark.optional
@pytest.mark.parametrize(
    "sample_key",
    ["sample", "sample_str"],
)
@pytest.mark.parametrize(
    "da_kwargs",
    [
        {"sample_cov_keys": ["meta1_cat"]},
        {"sample_cov_keys": ["meta1_cat", "batch"]},
        {"sample_cov_keys": ["meta1_cat"], "omit_original_sample": False},
        {"sample_cov_keys": ["meta1_cat"], "compute_log_enrichment": True},
        {"sample_cov_keys": ["meta1_cat", "batch"], "compute_log_enrichment": True},
    ],
)
def test_mrvi_da(model, sample_key, da_kwargs):
    model.differential_abundance(**da_kwargs)


@pytest.mark.optional
@pytest.mark.parametrize(
    "model_kwargs",
    [
        {"qz_kwargs": {"use_map": False}},
        {
            "qz_kwargs": {"use_map": False},
            "px_kwargs": {"low_dim_batch": False},
            "u_prior_mixture": True,
        },
        {
            "qz_kwargs": {
                "use_map": False,
                "stop_gradients": False,
                "stop_gradients_mlp": True,
            },
            "px_kwargs": {
                "low_dim_batch": False,
                "stop_gradients": False,
                "stop_gradients_mlp": True,
            },
            "z_u_prior": False,
        },
        {
            "qz_kwargs": {"use_map": False},
            "px_kwargs": {"low_dim_batch": True},
            "learn_z_u_prior_scale": True,
        },
    ],
)
def test_mrvi_model_kwargs(adata: AnnData, model_kwargs: dict[str, Any], save_path: str):
    MRVI.setup_anndata(
        adata,
        sample_key="sample_str",
        batch_key="batch",
    )
    model = MRVI(adata, n_latent=10, scale_observations=True, **model_kwargs)
    model.train(max_steps=2, train_size=0.5)

    model_path = os.path.join(save_path, "mrvi_model")
    model.save(model_path, save_anndata=False, overwrite=True)
    model = MRVI.load(model_path, adata=adata)


@pytest.mark.optional
def test_mrvi_sample_subset(model: MRVI, adata: AnnData, save_path: str):
    sample_cov_keys = ["meta1_cat", "meta2", "cont_cov"]
    sample_subset = [chr(i + ord("a")) for i in range(8)]
    model.differential_expression(sample_cov_keys=sample_cov_keys, sample_subset=sample_subset)

    model_path = os.path.join(save_path, "mrvi_model")
    model.save(model_path, save_anndata=False, overwrite=True)
    model = MRVI.load(model_path, adata=adata)


@pytest.mark.optional
def test_mrvi_shrink_u(adata: AnnData, save_path: str):
    MRVI.setup_anndata(
        adata,
        sample_key="sample_str",
        batch_key="batch",
    )
    model = MRVI(adata, n_latent=10, n_latent_u=5)
    model.train(max_steps=2, train_size=0.5)
    model.get_local_sample_distances()

    assert model.get_latent_representation().shape == (
        adata.shape[0],
        5,
    )

    model_path = os.path.join(save_path, "mrvi_model")
    model.save(model_path, save_anndata=False, overwrite=True)
    model = MRVI.load(model_path, adata=adata)


@pytest.fixture
def adata_stratifications():
    adata = synthetic_iid()
    adata.obs["sample"] = np.random.choice(15, size=adata.shape[0])
    adata.obs["sample_str"] = [chr(i + ord("a")) for i in adata.obs["sample"]]
    meta1 = np.random.randint(0, 2, size=15)
    adata.obs["meta1"] = meta1[adata.obs["sample"].values]
    meta2 = np.random.randn(15)
    adata.obs["meta2"] = meta2[adata.obs["sample"].values]
    adata.obs["cont_cov"] = np.random.normal(0, 1, size=adata.shape[0])
    adata.obs.loc[:, "label_2"] = np.random.choice(2, size=adata.shape[0])
    return adata


@pytest.mark.optional
def test_mrvi_stratifications(adata_stratifications: AnnData, save_path: str):
    MRVI.setup_anndata(
        adata_stratifications,
        sample_key="sample_str",
        batch_key="batch",
    )
    model = MRVI(adata_stratifications, n_latent=10)
    model.train(max_steps=2, train_size=0.5)

    dists = model.get_local_sample_distances(groupby=["labels", "label_2"])
    cell_dists = dists["cell"]
    assert cell_dists.shape == (adata_stratifications.shape[0], 15, 15)
    ct_dists = dists["labels"]
    assert ct_dists.shape == (3, 15, 15)
    assert np.allclose(ct_dists[0].values, ct_dists[0].values.T, atol=1e-6)
    ct_dists = dists["label_2"]
    assert ct_dists.shape == (2, 15, 15)
    assert np.allclose(ct_dists[0].values, ct_dists[0].values.T, atol=1e-6)

    model_path = os.path.join(save_path, "mrvi_model")
    model.save(model_path, save_anndata=False, overwrite=True)
    model = MRVI.load(model_path, adata=adata_stratifications)
