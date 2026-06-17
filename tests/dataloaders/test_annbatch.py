"""Annbatch disk-based dataloader tests for all supported scvi-tools models."""

from __future__ import annotations

import json
import os
import subprocess
import sys

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
    inference_dl = dm.inference_dataloader()

    de = model.differential_expression(
        groupby="labels",
        group1="label_1",
        pseudocounts=1e-4,
        n_samples_overall=50,
        silent=True,
        dataloader=inference_dl,
    )
    assert not de.empty
    assert {"bayes_factor", "group1", "group2"}.issubset(de.columns)

    de_change = model.differential_expression(
        groupby="labels",
        group1="label_1",
        pseudocounts=1e-4,
        n_samples_overall=50,
        silent=True,
        mode="change",
        dataloader=inference_dl,
    )
    assert not de_change.empty
    assert {"delta", "lfc_mean"}.issubset(de_change.columns)

    de_importance = model.differential_expression(
        groupby="labels",
        group1="label_1",
        weights="importance",
        n_samples_overall=50,
        silent=True,
        dataloader=inference_dl,
    )
    assert not de_importance.empty

    da = model.differential_abundance(
        sample_key="sample",
        batch_size=256,
        num_cells_posterior=50,
        dof=3,
        dataloader=inference_dl,
    )
    assert isinstance(da, pd.DataFrame)
    assert da.shape == (dm.n_obs, dm.n_samples)


@pytest.mark.dataloader
def test_annbatch_memory_contract(save_path: str):
    """setup_annbatch uses far less steady-state RSS than setup_anndata."""
    paths, _ = _synthetic_files(save_path, "annbatch_memory", batch_size=2000)
    collection_path = os.path.join(save_path, "annbatch_memory_collection.zarr")
    scvi.model.SCVI.setup_annbatch(
        collection_path=collection_path,
        paths=paths,
        batch_key="batch",
        labels_key="labels",
        batch_size=256,
    )

    script = """
import gc
import json
import sys

import anndata as ad
import psutil
import scvi

mode = sys.argv[1]
paths = json.loads(sys.argv[2])
collection_path = sys.argv[3]
proc = psutil.Process()


def rss():
    gc.collect()
    return proc.memory_info().rss


before = rss()
if mode == "anndata":
    adata = ad.concat([ad.read_h5ad(path) for path in paths], axis=0)
    scvi.model.SCVI.setup_anndata(adata, batch_key="batch", labels_key="labels")
    model = scvi.model.SCVI(adata)
elif mode == "annbatch":
    dm = scvi.model.SCVI.setup_annbatch(
        collection_path=collection_path,
        batch_key="batch",
        labels_key="labels",
        batch_size=256,
        rebuild=False,
    )
    model = scvi.model.SCVI(registry=dm.registry)
else:
    raise SystemExit(mode)

after = rss()
print(
    json.dumps(
        {
            "delta": after - before,
            "final": after,
            "has_adata": model.adata is not None,
        }
    )
)
"""

    def _probe(mode: str) -> dict[str, int | bool]:
        out = subprocess.check_output(
            [sys.executable, "-c", script, mode, json.dumps(paths), collection_path],
            text=True,
        )
        # Scan from the end for the JSON line; the subprocess may emit trailing
        # ANSI escape codes (e.g. a `\x1b[0m` reset) after the payload.
        for line in reversed(out.strip().splitlines()):
            line = line.strip()
            if line.startswith("{"):
                return json.loads(line)
        raise ValueError(f"No JSON payload in probe output:\n{out}")

    anndata_probe = _probe("anndata")
    annbatch_probe = _probe("annbatch")

    assert anndata_probe["has_adata"] is True
    assert annbatch_probe["has_adata"] is False
    assert annbatch_probe["delta"] < anndata_probe["delta"]
    assert anndata_probe["delta"] - annbatch_probe["delta"] > 10 * 1024 * 1024


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
    for mode in ("train", "validation"):
        for metric in ("classification_loss", "calibration_error", "accuracy"):
            assert f"{mode}_{metric}" in keys

    inference_dl = dm.inference_dataloader()
    latent = model.get_latent_representation(dataloader=inference_dl)
    assert latent.shape[0] == dm.n_obs

    predictions = model.predict(dataloader=inference_dl, soft=False)
    assert len(predictions) == dm.n_obs

    soft_preds = model.predict(dataloader=inference_dl, soft=True)
    assert soft_preds.shape[0] == dm.n_obs
    assert soft_preds.shape[1] == dm.n_labels

    _assert_save_load(model, scvi.model.SCANVI, save_path, "setup_scanvi", dm)


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


@pytest.mark.dataloader
@pytest.mark.parametrize("model_cls", [scvi.model.LinearSCVI, scvi.external.ContrastiveVI])
def test_annbatch_setup_dense_layer(
    model_cls,
    save_path: str,
):
    """setup_annbatch reads dense h5ad layers from the zarr collection."""
    _zarr()
    paths = []
    for i in range(2):
        adata = scvi.data.synthetic_iid(batch_size=20, n_genes=10)
        dense_counts = np.asarray(adata.X, dtype=np.float32)
        adata.X = dense_counts
        adata.layers["counts"] = dense_counts.copy()
        p = os.path.join(save_path, f"{model_cls.__name__.lower()}_dense_{i + 1}.h5ad")
        adata.write(p)
        paths.append(p)

    dm = model_cls.setup_annbatch(
        collection_path=os.path.join(save_path, f"{model_cls.__name__.lower()}_dense.zarr"),
        paths=paths,
        layer="counts",
        batch_size=16,
        chunk_size=8,
        preload_nchunks=2,
        dataset_size=1024,
    )

    assert dm.registry["setup_args"]["layer"] == "counts"
    batch = next(iter(dm.inference_dataloader()))
    assert batch["X"].shape == (16, dm.n_vars)


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

    ambient_profile = model.get_ambient_profile(
        dataloader=inference_dl,
        prob=0.0,
        iterations=1,
        sample=100,
    )
    assert ambient_profile.shape[0] == dm.n_vars
    assert ambient_profile.ndim == 2
    assert np.isfinite(ambient_profile).all()


@pytest.mark.dataloader
def test_annbatch_setup_mrvi(save_path: str):
    """TorchMRVI.setup_annbatch: build, train, latent, normalized_expression."""
    _zarr()
    paths = []
    donors_sites = [
        ("donor_A", "site_1"),
        ("donor_B", "site_1"),
        ("donor_C", "site_2"),
        ("donor_D", "site_2"),
    ]
    for i, (donor, site) in enumerate(donors_sites):
        adata = scvi.data.synthetic_iid(batch_size=250)
        adata.X = csr_matrix(adata.X)
        adata.obs["donor"] = donor
        adata.obs["site"] = pd.Categorical([site] * adata.n_obs, categories=["site_1", "site_2"])
        p = os.path.join(save_path, f"mrvi_{i + 1}.h5ad")
        adata.write(p)
        paths.append(p)

    collection_path = os.path.join(save_path, "mrvi.zarr")
    dm = scvi.external.TorchMRVI.setup_annbatch(
        collection_path=collection_path,
        paths=paths,
        batch_key="batch",
        sample_key="donor",
        categorical_covariate_keys=["site"],
        batch_size=256,
        dataset_size=1024,
    )
    assert dm.n_batch == 2
    assert dm.n_samples == 4
    assert dm.registry["field_registries"]["sample"]["summary_stats"]["n_sample"] == 4
    assert (
        dm.registry["field_registries"]["extra_categorical_covs"]["summary_stats"][
            "n_extra_categorical_covs"
        ]
        == 1
    )

    model = scvi.external.TorchMRVI(registry=dm.registry)
    model.train(max_epochs=1, datamodule=dm, train_size=0.8, check_val_every_n_epoch=1)
    assert "elbo_train" in model.history
    _assert_validation_split(model, dm)

    inference_dl = dm.inference_dataloader()

    latent = model.get_latent_representation(dataloader=inference_dl, give_z=False)
    assert latent.shape[0] == dm.n_obs

    norm_expr = model.get_normalized_expression(dataloader=inference_dl)
    assert norm_expr.shape[0] == dm.n_obs

    model.update_sample_info(dataloader=inference_dl)
    assert "site" in model.sample_info.columns
    assert not model.sample_info.empty

    local_sample_distances = model.get_local_sample_distances(
        dataloader=inference_dl,
        batch_size=128,
    )
    assert "cell" in local_sample_distances.data_vars


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
        labels_key="labels",
        batch_size=256,
        dataset_size=1024,
    )
    assert dm.n_batch == 2
    assert dm.n_labels > 1

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

    da = model.differential_accessibility(
        groupby="labels",
        group1=dm.labels_[1],
        batch_size=256,
        dataloader=inference_dl,
        silent=True,
    )
    assert not da.empty
    assert "bayes_factor" in da.columns

    reconstruction = model.get_reconstruction_error(dataloader=inference_dl)
    assert isinstance(reconstruction, dict)

    # Region factors — per-region parameter stored in module, no adata needed
    region_factors = model.get_region_factors()
    assert len(region_factors) == dm.n_vars

    _assert_save_load(model, scvi.model.PEAKVI, save_path, "setup_peakvi", dm)


@pytest.mark.dataloader
def test_annbatch_in_memory_basic(save_path: str):
    """In-memory adatas path trains SCVI without writing zarr."""
    adatas = []
    for _i in range(2):
        adata = scvi.data.synthetic_iid(batch_size=500, n_genes=20)
        adata.X = csr_matrix(adata.X)
        adatas.append(adata)

    dm = scvi.model.SCVI.setup_annbatch(
        adatas=adatas,
        batch_key="batch",
        labels_key="labels",
        batch_size=256,
        chunk_size=64,
        preload_nchunks=4,
    )

    assert dm.n_batch == 2
    assert dm.n_vars == adatas[0].n_vars
    col_names = dm.registry["field_registries"]["X"]["state_registry"]["column_names"]
    assert col_names == list(adatas[0].var_names)

    model = scvi.model.SCVI(registry=dm.registry)
    model.train(max_epochs=1, datamodule=dm, train_size=0.8, check_val_every_n_epoch=1)
    assert "elbo_train" in model.history

    inference_dl = dm.inference_dataloader()
    latent = model.get_latent_representation(dataloader=inference_dl)
    assert latent.shape == (dm.n_obs, model.module.n_latent)

    _assert_save_load(model, scvi.model.SCVI, save_path, "in_memory_basic", dm)


@pytest.mark.dataloader
def test_annbatch_in_memory_collection_path_warns(save_path: str):
    """Passing both adatas and collection_path warns and ignores collection_path."""
    adatas = [scvi.data.synthetic_iid(batch_size=50, n_genes=5) for _ in range(2)]
    for a in adatas:
        a.X = csr_matrix(a.X)

    with pytest.warns(UserWarning, match="collection_path.*ignored"):
        dm = scvi.model.SCVI.setup_annbatch(
            adatas=adatas,
            collection_path=os.path.join(save_path, "should_not_exist.zarr"),
            batch_key="batch",
            batch_size=16,
            chunk_size=8,
            preload_nchunks=2,
        )
    assert not os.path.exists(os.path.join(save_path, "should_not_exist.zarr"))
    assert dm.n_obs == sum(a.n_obs for a in adatas)


@pytest.mark.dataloader
def test_annbatch_in_memory_paths_error():
    """Passing both adatas and paths raises ValueError."""
    adatas = [scvi.data.synthetic_iid(batch_size=10, n_genes=5)]
    with pytest.raises(ValueError, match="mutually exclusive"):
        scvi.model.SCVI.setup_annbatch(
            adatas=adatas,
            paths=["some_file.h5ad"],
            batch_key="batch",
        )


@pytest.mark.dataloader
def test_annbatch_class_sampler_labels(save_path: str):
    """ClassSampler is used for training when use_class_sampler=True with labels_key."""
    _zarr()
    # Each file contains a single label → contiguous runs satisfy ClassSampler's run-length rule
    paths = []
    for i, label in enumerate(["label_0", "label_1"]):
        adata = scvi.data.synthetic_iid(batch_size=500, n_genes=10)
        adata.X = csr_matrix(adata.X)
        adata.obs["labels"] = label
        p = os.path.join(save_path, f"cs_labels_{i + 1}.h5ad")
        adata.write(p)
        paths.append(p)

    dm = scvi.model.SCVI.setup_annbatch(
        collection_path=os.path.join(save_path, "cs_labels.zarr"),
        paths=paths,
        batch_key="batch",
        labels_key="labels",
        use_class_sampler=True,
        batch_size=64,
        chunk_size=64,
        preload_nchunks=4,
        dataset_size=1024,
    )

    assert dm._classes is not None
    assert set(dm._classes.categories) == {"label_0", "label_1"}

    # After set_split, _shuffle_train=True so train_dataloader uses ClassSampler
    from annbatch.samplers import ClassSampler

    dm.set_split(train_size=0.9)
    train_loader = dm.train_dataloader()
    assert isinstance(train_loader._batch_sampler, ClassSampler)

    # Verify each batch is pure-class
    for raw_batch in train_loader:
        obs = raw_batch.get("obs")
        assert obs is not None
        assert obs["labels"].nunique() == 1, "Expected pure-class batch from ClassSampler"
        break

    model = scvi.model.SCVI(registry=dm.registry)
    model.train(max_epochs=1, datamodule=dm, train_size=0.9, check_val_every_n_epoch=1)
    assert "elbo_train" in model.history


@pytest.mark.dataloader
def test_annbatch_class_sampler_custom_key(save_path: str):
    """class_sampler_key overrides labels_key for class assignment."""
    _zarr()
    paths = []
    for i, batch_val in enumerate(["batch_0", "batch_1"]):
        adata = scvi.data.synthetic_iid(batch_size=200, n_genes=10)
        adata.X = csr_matrix(adata.X)
        adata.obs["batch"] = batch_val
        p = os.path.join(save_path, f"cs_custom_{i + 1}.h5ad")
        adata.write(p)
        paths.append(p)

    dm = scvi.model.SCVI.setup_annbatch(
        collection_path=os.path.join(save_path, "cs_custom.zarr"),
        paths=paths,
        batch_key="batch",
        use_class_sampler=True,
        class_sampler_key="batch",
        batch_size=64,
        chunk_size=64,
        preload_nchunks=4,
        dataset_size=1024,
    )

    assert dm._classes is not None
    assert set(dm._classes.categories) == {"batch_0", "batch_1"}

    from annbatch.samplers import ClassSampler

    dm.set_split(train_size=0.9)
    assert isinstance(dm.train_dataloader()._batch_sampler, ClassSampler)


@pytest.mark.dataloader
def test_annbatch_class_sampler_no_key_error():
    """use_class_sampler=True with no labels_key or class_sampler_key raises ValueError."""
    with pytest.raises(ValueError, match="requires.*labels_key.*class_sampler_key"):
        scvi.model.SCVI.setup_annbatch(
            paths=["dummy.h5ad"],
            use_class_sampler=True,
            batch_size=64,
        )


@pytest.mark.dataloader
def test_annbatch_class_sampler_shuffle_warns(save_path: str):
    """use_class_sampler=True with shuffle=True emits a UserWarning."""
    import contextlib

    paths = []
    for i, label in enumerate(["label_0", "label_1"]):
        adata = scvi.data.synthetic_iid(batch_size=50, n_genes=5)
        adata.X = csr_matrix(adata.X)
        adata.obs["labels"] = label
        p = os.path.join(save_path, f"cs_shuffle_warn_{i + 1}.h5ad")
        adata.write(p)
        paths.append(p)

    with pytest.warns(UserWarning, match="shuffle.*class locality"):
        with contextlib.suppress(Exception):
            scvi.model.SCVI.setup_annbatch(
                collection_path=os.path.join(save_path, "cs_shuffle_warn.zarr"),
                paths=paths,
                labels_key="labels",
                use_class_sampler=True,
                shuffle=True,
                batch_size=16,
                chunk_size=8,
                preload_nchunks=2,
                dataset_size=256,
            )


@pytest.mark.dataloader
def test_annbatch_in_memory_class_sampler_error():
    """use_class_sampler=True with adatas raises ValueError (annbatch limitation)."""
    adatas = [scvi.data.synthetic_iid(batch_size=50, n_genes=5) for _ in range(2)]
    for a in adatas:
        a.X = csr_matrix(a.X)

    with pytest.raises(ValueError, match="ClassSampler.*not supported.*adatas"):
        scvi.model.SCVI.setup_annbatch(
            adatas=adatas,
            labels_key="labels",
            use_class_sampler=True,
            batch_size=16,
        )
