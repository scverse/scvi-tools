"""Annbatch disk-based dataloader tests for all supported scvi-tools models."""

from __future__ import annotations

import os

import numpy as np
import pandas as pd
import pytest
from scipy.sparse import csr_matrix

import scvi

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------


def _zarr():
    import zarr

    zarr.config.set({"codec_pipeline.path": "zarrs.ZarrsCodecPipeline"})


def _synthetic_files(save_path, prefix, n=2, batch_size=500, **obs_cols):
    """Write n synthetic h5ad files; return (paths, last_adata)."""
    paths = []
    last = None
    for i in range(n):
        adata = scvi.data.synthetic_iid(batch_size=batch_size)
        adata.X = csr_matrix(adata.X)
        for col, vals in obs_cols.items():
            adata.obs[col] = vals[i] if isinstance(vals, list) else vals
        p = os.path.join(save_path, f"{prefix}_{i + 1}.h5ad")
        adata.write(p)
        paths.append(p)
        last = adata
    return paths, last


def _assert_validation_split(model, dm):
    """Validate that train_size created a non-overlapping validation stream."""
    assert dm.n_train > 0
    assert dm.n_val > 0
    assert len(dm.train_idx) == dm.n_train
    assert len(dm.val_idx) == dm.n_val
    assert np.intersect1d(dm.train_idx, dm.val_idx).size == 0
    assert any("validation" in key for key in model.history.keys())


def _assert_save_load(model, model_cls, save_path, model_name, dm):
    """Save annbatch-trained model and reload it, verifying inference works."""
    save_dir = os.path.join(save_path, f"saved_{model_name}")
    # save_anndata=True is a no-op for annbatch models (self.adata is None)
    model.save(save_dir, overwrite=True, datamodule=dm)

    model_loaded = model_cls.load(save_dir, adata=False)
    assert model_loaded.is_trained_

    if hasattr(model_loaded, "get_latent_representation"):
        inference_dl = dm.inference_dataloader()
        latent = model_loaded.get_latent_representation(dataloader=inference_dl)
        assert latent.shape[0] == dm.n_obs


# ---------------------------------------------------------------------------
# Generic annbatch / SCVI
# ---------------------------------------------------------------------------


@pytest.mark.dataloader
def test_annbatch(save_path: str):
    """Basic SCVI annbatch: build, train, infer, validate."""
    paths, ref = _synthetic_files(save_path, "annbatch_basic", batch_size=20000)

    dm = scvi.model.SCVI.setup_annbatch(
        paths=paths,
        batch_key="batch",
        labels_key="labels",
    )

    assert dm.n_batch == 2
    assert dm.n_labels > 0

    model = scvi.model.SCVI(registry=dm.registry)
    model.train(max_epochs=1, datamodule=dm, train_size=0.8, check_val_every_n_epoch=1)
    assert "elbo_train" in model.history
    _assert_validation_split(model, dm)

    inference_dl = dm.inference_dataloader()
    latent = model.get_latent_representation(dataloader=inference_dl)
    assert latent.shape == (dm.n_obs, model.module.n_latent)
    _ = model.get_elbo(dataloader=inference_dl)

    _assert_save_load(model, scvi.model.SCVI, save_path, "annbatch_basic", dm)


@pytest.mark.dataloader
def test_annbatch_with_covariates(save_path: str):
    """SCVI annbatch with categorical and continuous covariates."""
    _zarr()
    paths, _ = _synthetic_files(save_path, "annbatch_cov", batch_size=20000)
    # Add covariates post-write by re-reading isn't needed — just write with obs cols
    # Rewrite with covariate columns
    paths = []
    for i in range(2):
        adata = scvi.data.synthetic_iid(batch_size=20000)
        adata.X = csr_matrix(adata.X)
        adata.obs["cat1"] = np.random.randint(0, 5, size=adata.n_obs).astype(str)
        adata.obs["cat2"] = np.random.randint(0, 3, size=adata.n_obs).astype(str)
        adata.obs["cont1"] = np.random.normal(size=adata.n_obs)
        adata.obs["cont2"] = np.random.normal(size=adata.n_obs)
        p = os.path.join(save_path, f"annbatch_cov_{i + 1}.h5ad")
        adata.write(p)
        paths.append(p)

    collection_path = os.path.join(save_path, "annbatch_covariates_collection")
    dm = scvi.model.SCVI.setup_annbatch(
        collection_path=collection_path,
        paths=paths,
        batch_key="batch",
        labels_key="labels",
        categorical_covariate_keys=["cat1", "cat2"],
        continuous_covariate_keys=["cont1", "cont2"],
        batch_size=4096,
        dataset_size=2_097_152,
    )

    reg = dm.registry
    assert (
        reg["field_registries"]["extra_categorical_covs"]["summary_stats"][
            "n_extra_categorical_covs"
        ]
        == 2
    )
    assert (
        reg["field_registries"]["extra_continuous_covs"]["summary_stats"][
            "n_extra_continuous_covs"
        ]
        == 2
    )

    model = scvi.model.SCVI(registry=dm.registry)
    model.train(max_epochs=1, datamodule=dm, train_size=0.8, check_val_every_n_epoch=1)
    _assert_validation_split(model, dm)

    keys = model.history.keys()
    assert "elbo_train" in keys
    assert "reconstruction_loss_train" in keys
    assert "kl_local_train" in keys

    inference_dl = dm.inference_dataloader()
    _ = model.get_elbo(dataloader=inference_dl)
    _ = model.get_latent_representation(dataloader=inference_dl)

    _assert_save_load(model, scvi.model.SCVI, save_path, "annbatch_cov", dm)


@pytest.mark.dataloader
def test_annbatch_scvi_downstream_tasks(save_path: str):
    """SCVI annbatch supports DE and differential abundance after training."""
    _zarr()
    paths = []
    for i in range(2):
        adata = scvi.data.synthetic_iid(batch_size=500, n_genes=20)
        adata.X = csr_matrix(adata.X)
        adata.obs["sample"] = np.where(np.arange(adata.n_obs) % 2 == 0, "sample_0", "sample_1")
        p = os.path.join(save_path, f"annbatch_downstream_{i + 1}.h5ad")
        adata.write(p)
        paths.append(p)

    dm = scvi.model.SCVI.setup_annbatch(
        collection_path=os.path.join(save_path, "annbatch_downstream_collection"),
        paths=paths,
        batch_key="batch",
        labels_key="labels",
        sample_key="sample",
        batch_size=256,
        dataset_size=1024,
    )

    model = scvi.model.SCVI(registry=dm.registry)
    model.train(max_epochs=1, datamodule=dm, train_size=0.8, check_val_every_n_epoch=1)
    _assert_validation_split(model, dm)

    de = model.differential_expression(
        groupby="labels",
        group1="label_1",
        pseudocounts=1e-4,
        n_samples_overall=50,
        silent=True,
    )
    assert not de.empty
    assert {"bayes_factor", "group1", "group2"}.issubset(de.columns)

    de_importance = model.differential_expression(
        groupby="labels",
        group1="label_1",
        weights="importance",
        n_samples_overall=50,
        silent=True,
    )
    assert not de_importance.empty

    model.differential_abundance(
        sample_key="sample",
        batch_size=256,
        num_cells_posterior=50,
        dof=3,
    )
    da = model.adata.obsm["da_log_probs"]
    assert isinstance(da, pd.DataFrame)
    assert da.shape == (dm.n_obs, model.adata.obs["sample"].nunique())


# ---------------------------------------------------------------------------
# SCVI setup_annbatch
# ---------------------------------------------------------------------------


@pytest.mark.dataloader
def test_annbatch_setup_scvi(save_path: str):
    """SCVI.setup_annbatch: build-once, reuse, rebuild, layer, validation."""
    _zarr()
    paths, ref = _synthetic_files(save_path, "scvi_setup", batch_size=500)
    collection_path = os.path.join(save_path, "annbatch_setup_scvi.zarr")

    # First build
    dm = scvi.model.SCVI.setup_annbatch(
        collection_path=collection_path,
        paths=paths,
        batch_key="batch",
        labels_key="labels",
        batch_size=256,
        dataset_size=1024,
    )

    assert dm.n_batch == 2
    assert dm.n_vars == ref.n_vars
    col_names = dm.registry["field_registries"]["X"]["state_registry"]["column_names"]
    assert col_names == list(ref.var_names)

    model = scvi.model.SCVI(registry=dm.registry)
    model.train(max_epochs=1, datamodule=dm, train_size=0.8, check_val_every_n_epoch=1)
    assert "elbo_train" in model.history
    _assert_validation_split(model, dm)

    inference_dl = dm.inference_dataloader()
    latent = model.get_latent_representation(dataloader=inference_dl)
    assert latent.shape == (dm.n_obs, model.module.n_latent)
    elbo = model.get_elbo(dataloader=inference_dl)
    assert isinstance(float(elbo), float)
    reconstruction = model.get_reconstruction_error(dataloader=inference_dl)
    assert isinstance(reconstruction, dict)

    # Reuse existing collection (no paths)
    dm2 = scvi.model.SCVI.setup_annbatch(
        collection_path=collection_path,
        batch_key="batch",
        batch_size=256,
        rebuild=False,
    )
    assert dm2.n_batch == 2
    assert dm2.n_vars == ref.n_vars

    # Explicit rebuild
    dm3 = scvi.model.SCVI.setup_annbatch(
        collection_path=collection_path,
        paths=paths,
        batch_key="batch",
        batch_size=256,
        rebuild=True,
    )
    model3 = scvi.model.SCVI(registry=dm3.registry)
    model3.train(max_epochs=1, datamodule=dm3)
    assert "elbo_train" in model3.history

    # With layer
    paths_layer = []
    for i, _base_path in enumerate(paths):
        adata_l = scvi.data.synthetic_iid(batch_size=500)
        adata_l.X = csr_matrix(adata_l.X)
        adata_l.layers["counts"] = adata_l.X.copy()
        p = os.path.join(save_path, f"scvi_setup_layer_{i + 1}.h5ad")
        adata_l.write(p)
        paths_layer.append(p)

    dm_layer = scvi.model.SCVI.setup_annbatch(
        collection_path=os.path.join(save_path, "annbatch_setup_scvi_layer.zarr"),
        paths=paths_layer,
        batch_key="batch",
        layer="counts",
        batch_size=256,
        dataset_size=1024,
    )
    model_layer = scvi.model.SCVI(registry=dm_layer.registry)
    model_layer.train(max_epochs=1, datamodule=dm_layer)
    assert "elbo_train" in model_layer.history

    _assert_save_load(model, scvi.model.SCVI, save_path, "setup_scvi", dm)


# ---------------------------------------------------------------------------
# SCANVI
# ---------------------------------------------------------------------------


@pytest.mark.dataloader
def test_annbatch_setup_scanvi(save_path: str):
    """SCANVI.setup_annbatch: train, latent, predict."""
    _zarr()
    paths, ref = _synthetic_files(save_path, "scanvi", batch_size=500)
    collection_path = os.path.join(save_path, "annbatch_scanvi.zarr")

    dm = scvi.model.SCANVI.setup_annbatch(
        collection_path=collection_path,
        paths=paths,
        batch_key="batch",
        labels_key="labels",
        unlabeled_category="Unknown",
        batch_size=256,
        dataset_size=1024,
    )

    assert dm.n_batch == 2
    assert dm.n_labels > 0
    assert dm.registry["model_name"] == "SCANVI"
    assert (
        dm.registry["field_registries"]["labels"]["state_registry"]["unlabeled_category"]
        == "Unknown"
    )

    model = scvi.model.SCANVI(
        adata=None, registry=dm.registry, encode_covariates=False, datamodule=dm
    )
    model.train(max_epochs=1, datamodule=dm, train_size=0.8, check_val_every_n_epoch=1)
    _assert_validation_split(model, dm)

    keys = model.history.keys()
    assert "elbo_train" in keys
    assert "train_classification_loss" in keys
    assert "train_accuracy" in keys

    inference_dl = dm.inference_dataloader()
    latent = model.get_latent_representation(dataloader=inference_dl)
    assert latent.shape[0] == dm.n_obs

    predictions = model.predict(dataloader=inference_dl, soft=False)
    assert len(predictions) == dm.n_obs

    soft_preds = model.predict(dataloader=inference_dl, soft=True)
    assert soft_preds.shape[0] == dm.n_obs
    assert soft_preds.shape[1] == dm.n_labels

    _assert_save_load(model, scvi.model.SCANVI, save_path, "setup_scanvi", dm)


# ---------------------------------------------------------------------------
# sample_key forwarding
# ---------------------------------------------------------------------------


@pytest.mark.dataloader
def test_annbatch_setup_base_sample_key(save_path: str):
    """Base setup_annbatch must forward sample_key to AnnbatchDataModule."""
    _zarr()
    paths = []
    for i, sample_name in enumerate(["sample_A", "sample_B"]):
        adata = scvi.data.synthetic_iid(batch_size=500)
        adata.X = csr_matrix(adata.X)
        adata.obs["sample"] = sample_name
        p = os.path.join(save_path, f"base_sample_{i + 1}.h5ad")
        adata.write(p)
        paths.append(p)

    dm = scvi.model.SCVI.setup_annbatch(
        collection_path=os.path.join(save_path, "base_sample.zarr"),
        paths=paths,
        batch_key="batch",
        sample_key="sample",
        batch_size=256,
        dataset_size=1024,
    )
    assert dm.n_samples == 2
    assert dm.registry["field_registries"]["sample"]["summary_stats"]["n_sample"] == 2


# ---------------------------------------------------------------------------
# LinearSCVI
# ---------------------------------------------------------------------------


@pytest.mark.dataloader
def test_annbatch_setup_linear_scvi(save_path: str):
    """LinearSCVI.setup_annbatch: build, train, latent, loadings."""
    _zarr()
    paths, ref = _synthetic_files(save_path, "linear_scvi", batch_size=500)
    collection_path = os.path.join(save_path, "linear_scvi.zarr")

    dm = scvi.model.LinearSCVI.setup_annbatch(
        collection_path=collection_path,
        paths=paths,
        batch_key="batch",
        batch_size=256,
        dataset_size=1024,
    )
    assert dm.n_batch == 2
    assert dm.n_vars == ref.n_vars

    model = scvi.model.LinearSCVI(registry=dm.registry)
    model.train(max_epochs=1, datamodule=dm, train_size=0.8, check_val_every_n_epoch=1)
    assert "elbo_train" in model.history
    _assert_validation_split(model, dm)

    inference_dl = dm.inference_dataloader()
    latent = model.get_latent_representation(dataloader=inference_dl)
    assert latent.shape[0] == dm.n_obs
    assert latent.shape[1] == model.n_latent

    # Loadings: decoder weight matrix — no adata needed
    loadings = model.get_loadings()
    assert loadings.shape == (dm.n_vars, model.n_latent)

    reconstruction = model.get_reconstruction_error(dataloader=inference_dl)
    assert isinstance(reconstruction, dict)

    _assert_save_load(model, scvi.model.LinearSCVI, save_path, "setup_linear_scvi", dm)


# ---------------------------------------------------------------------------
# AUTOZI
# ---------------------------------------------------------------------------


@pytest.mark.dataloader
def test_annbatch_setup_autozi(save_path: str):
    """AUTOZI.setup_annbatch: build, train, elbo, reconstruction, alphas_betas."""
    _zarr()
    paths, _ = _synthetic_files(save_path, "autozi", batch_size=500)
    collection_path = os.path.join(save_path, "autozi.zarr")

    dm = scvi.model.AUTOZI.setup_annbatch(
        collection_path=collection_path,
        paths=paths,
        batch_key="batch",
        batch_size=256,
        dataset_size=1024,
    )
    assert dm.n_batch == 2

    model = scvi.model.AUTOZI(registry=dm.registry)
    model.train(max_epochs=1, datamodule=dm, train_size=0.8, check_val_every_n_epoch=1)
    assert "elbo_train" in model.history
    _assert_validation_split(model, dm)

    inference_dl = dm.inference_dataloader()
    elbo = model.get_elbo(dataloader=inference_dl)
    assert isinstance(float(elbo), float)

    reconstruction = model.get_reconstruction_error(dataloader=inference_dl)
    assert isinstance(reconstruction, dict)

    # Model-wide learned parameters — no adata needed
    alphas_betas = model.get_alphas_betas()
    assert "alpha_posterior" in alphas_betas
    assert "beta_posterior" in alphas_betas

    _assert_save_load(model, scvi.model.AUTOZI, save_path, "setup_autozi", dm)


# ---------------------------------------------------------------------------
# CondSCVI
# ---------------------------------------------------------------------------


@pytest.mark.dataloader
def test_annbatch_setup_condscvi(save_path: str):
    """CondSCVI.setup_annbatch: build, train, latent, elbo, reconstruction."""
    _zarr()
    paths, _ = _synthetic_files(save_path, "condscvi", batch_size=500)
    collection_path = os.path.join(save_path, "condscvi.zarr")

    dm = scvi.model.CondSCVI.setup_annbatch(
        collection_path=collection_path,
        paths=paths,
        batch_key="batch",
        labels_key="labels",
        batch_size=256,
        dataset_size=1024,
    )
    assert dm.n_batch == 2
    assert dm.n_labels > 0

    model = scvi.model.CondSCVI(registry=dm.registry)
    model.train(max_epochs=1, datamodule=dm, train_size=0.8, check_val_every_n_epoch=1)
    assert "elbo_train" in model.history
    _assert_validation_split(model, dm)

    inference_dl = dm.inference_dataloader()
    latent = model.get_latent_representation(dataloader=inference_dl)
    assert latent.shape[0] == dm.n_obs

    elbo = model.get_elbo(dataloader=inference_dl)
    assert isinstance(float(elbo), float)

    reconstruction = model.get_reconstruction_error(dataloader=inference_dl)
    assert isinstance(reconstruction, dict)

    _assert_save_load(model, scvi.model.CondSCVI, save_path, "setup_condscvi", dm)


# ---------------------------------------------------------------------------
# SysVI
# ---------------------------------------------------------------------------


@pytest.mark.dataloader
def test_annbatch_setup_sysvi(save_path: str):
    """SysVI.setup_annbatch: build, train, latent, normalized_expression."""
    _zarr()
    paths, _ = _synthetic_files(save_path, "sysvi", batch_size=500)
    collection_path = os.path.join(save_path, "sysvi.zarr")

    dm = scvi.external.SysVI.setup_annbatch(
        collection_path=collection_path,
        paths=paths,
        batch_key="batch",
        batch_size=256,
        dataset_size=1024,
    )
    assert dm.n_batch == 2

    model = scvi.external.SysVI(registry=dm.registry, prior="standard_normal")
    model.train(max_epochs=1, datamodule=dm, train_size=0.8, check_val_every_n_epoch=1)
    assert "elbo_train" in model.history
    _assert_validation_split(model, dm)

    inference_dl = dm.inference_dataloader()
    latent = model.get_latent_representation(dataloader=inference_dl)
    assert latent.shape[0] == dm.n_obs

    norm_expr = model.get_normalized_expression(dataloader=inference_dl, library_size="latent")
    assert norm_expr.shape == (dm.n_obs, dm.n_vars)

    _assert_save_load(model, scvi.external.SysVI, save_path, "setup_sysvi", dm)


# ---------------------------------------------------------------------------
# SCAR
# ---------------------------------------------------------------------------


@pytest.mark.dataloader
def test_annbatch_setup_scar(save_path: str):
    """SCAR.setup_annbatch: build, train, elbo, latent, reconstruction, marginal_ll."""
    _zarr()
    paths, ref = _synthetic_files(save_path, "scar", batch_size=500)
    collection_path = os.path.join(save_path, "scar.zarr")

    dm = scvi.external.SCAR.setup_annbatch(
        collection_path=collection_path,
        paths=paths,
        batch_key="batch",
        batch_size=256,
        dataset_size=1024,
    )
    assert dm.n_batch == 2

    n_genes = ref.n_vars
    ambient_profile = np.full((1, n_genes), 1.0 / n_genes, dtype=np.float32)
    model = scvi.external.SCAR(registry=dm.registry, ambient_profile=ambient_profile)
    model.train(max_epochs=1, datamodule=dm, train_size=0.8, check_val_every_n_epoch=1)
    assert "elbo_train" in model.history
    _assert_validation_split(model, dm)

    inference_dl = dm.inference_dataloader()

    elbo = model.get_elbo(dataloader=inference_dl)
    assert isinstance(float(elbo), float)

    latent = model.get_latent_representation(dataloader=inference_dl)
    assert latent.shape[0] == dm.n_obs

    reconstruction = model.get_reconstruction_error(dataloader=inference_dl)
    assert isinstance(reconstruction, dict)

    marginal_ll = model.get_marginal_ll(dataloader=inference_dl, n_mc_samples=3)
    assert isinstance(marginal_ll, float)

    _assert_save_load(model, scvi.external.SCAR, save_path, "setup_scar", dm)


# ---------------------------------------------------------------------------
# MRVI / TorchMRVI
# ---------------------------------------------------------------------------


@pytest.mark.dataloader
def test_annbatch_setup_mrvi(save_path: str):
    """TorchMRVI.setup_annbatch: build, train, latent, normalized_expression."""
    _zarr()
    paths = []
    for i, donor in enumerate(["donor_A", "donor_B"]):
        adata = scvi.data.synthetic_iid(batch_size=500)
        adata.X = csr_matrix(adata.X)
        adata.obs["donor"] = donor
        p = os.path.join(save_path, f"mrvi_{i + 1}.h5ad")
        adata.write(p)
        paths.append(p)

    collection_path = os.path.join(save_path, "mrvi.zarr")
    dm = scvi.external.TorchMRVI.setup_annbatch(
        collection_path=collection_path,
        paths=paths,
        batch_key="batch",
        sample_key="donor",
        batch_size=256,
        dataset_size=1024,
    )
    assert dm.n_batch == 2
    assert dm.n_samples == 2
    assert dm.registry["field_registries"]["sample"]["summary_stats"]["n_sample"] == 2

    model = scvi.external.TorchMRVI(registry=dm.registry)
    model.train(max_epochs=1, datamodule=dm, train_size=0.8, check_val_every_n_epoch=1)
    assert "elbo_train" in model.history
    _assert_validation_split(model, dm)

    inference_dl = dm.inference_dataloader()

    latent = model.get_latent_representation(dataloader=inference_dl, give_z=False)
    assert latent.shape[0] == dm.n_obs

    norm_expr = model.get_normalized_expression(dataloader=inference_dl)
    assert norm_expr.shape[0] == dm.n_obs

    _assert_save_load(model, scvi.external.TorchMRVI, save_path, "setup_mrvi", dm)


# ---------------------------------------------------------------------------
# PEAKVI
# ---------------------------------------------------------------------------


@pytest.mark.dataloader
def test_annbatch_setup_peakvi(save_path: str):
    """PEAKVI.setup_annbatch: build, train, latent, accessibility, region_factors."""
    _zarr()
    paths = []
    for i in range(2):
        adata = scvi.data.synthetic_iid(batch_size=500)
        adata.X = csr_matrix((adata.X > 0).astype(np.float32))
        p = os.path.join(save_path, f"peakvi_{i + 1}.h5ad")
        adata.write(p)
        paths.append(p)

    collection_path = os.path.join(save_path, "peakvi.zarr")
    dm = scvi.model.PEAKVI.setup_annbatch(
        collection_path=collection_path,
        paths=paths,
        batch_key="batch",
        batch_size=256,
        dataset_size=1024,
    )
    assert dm.n_batch == 2

    model = scvi.model.PEAKVI(registry=dm.registry)
    model.train(
        max_epochs=1,
        datamodule=dm,
        train_size=0.8,
        check_val_every_n_epoch=1,
        early_stopping=False,
    )
    assert "elbo_train" in model.history
    _assert_validation_split(model, dm)

    inference_dl = dm.inference_dataloader()

    latent = model.get_latent_representation(dataloader=inference_dl)
    assert latent.shape[0] == dm.n_obs

    accessibility = model.get_normalized_accessibility(dataloader=inference_dl)
    assert accessibility.shape == (dm.n_obs, dm.n_vars)

    reconstruction = model.get_reconstruction_error(dataloader=inference_dl)
    assert isinstance(reconstruction, dict)

    # Region factors — per-region parameter stored in module, no adata needed
    region_factors = model.get_region_factors()
    assert len(region_factors) == dm.n_vars

    _assert_save_load(model, scvi.model.PEAKVI, save_path, "setup_peakvi", dm)
