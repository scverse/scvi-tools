"""Comprehensive tests for the ScPoli model."""

from __future__ import annotations

import os

import numpy as np
import pandas as pd
import pytest
import torch

import scvi
from scvi.data import synthetic_iid
from scvi.data._constants import ADATA_MINIFY_TYPE
from scvi.data._utils import _is_minified
from scvi.external import ScPoli

# ──────────────────────────────────────────────────────────────────────────────
# Constants / shared helpers
# ──────────────────────────────────────────────────────────────────────────────

_UNLABELED = "label_0"
_SCPOLI_OBSERVED_LIB_SIZE = "_scpoli_observed_lib_size"


def _setup_and_train(
    adata=None,
    labels_key: str = "labels",
    unlabeled_category: str = _UNLABELED,
    batch_key: str = "batch",
    layer: str | None = None,
    size_factor_key: str | None = None,
    categorical_covariate_keys: list[str] | None = None,
    continuous_covariate_keys: list[str] | None = None,
    n_latent: int = 5,
    max_epochs: int = 3,
    pretrain_epochs: int = 2,
    train_size: float = 0.9,
    **model_kwargs,
) -> tuple[ScPoli, object]:
    """Build, register and train a minimal ScPoli model on synthetic data."""
    if adata is None:
        adata = synthetic_iid()
    ScPoli.setup_anndata(
        adata,
        labels_key=labels_key,
        unlabeled_category=unlabeled_category,
        batch_key=batch_key,
        layer=layer,
        size_factor_key=size_factor_key,
        categorical_covariate_keys=categorical_covariate_keys,
        continuous_covariate_keys=continuous_covariate_keys,
    )
    model = ScPoli(adata, n_latent=n_latent, **model_kwargs)
    model.train(
        max_epochs=max_epochs,
        pretrain_epochs=pretrain_epochs,
        train_size=train_size,
    )
    return model, adata


# ──────────────────────────────────────────────────────────────────────────────
# Core functionality
# ──────────────────────────────────────────────────────────────────────────────


def test_scpoli():
    """Basic smoke test: indices, training phases, history keys, predict, prototypes."""
    adata = synthetic_iid()
    ScPoli.setup_anndata(
        adata, labels_key="labels", unlabeled_category=_UNLABELED, batch_key="batch"
    )
    model = ScPoli(adata, n_latent=5)

    # Labeled / unlabeled split is populated at __init__
    assert len(model._labeled_indices) == sum(adata.obs["labels"] != _UNLABELED)
    assert len(model._unlabeled_indices) == sum(adata.obs["labels"] == _UNLABELED)
    assert model.n_labels == len(adata.obs["labels"].unique()) - 1  # minus unlabeled

    # Prototypes should not yet be initialized
    assert not model.module._prototypes_initialized.item()

    model.train(max_epochs=3, pretrain_epochs=2, train_size=0.9, check_val_every_n_epoch=1)

    # After training past pretrain_epochs, prototypes are initialized
    assert model.module._prototypes_initialized.item()

    # History keys
    logged = model.history.keys()
    assert "elbo_train" in logged
    assert "reconstruction_loss_train" in logged
    assert "kl_local_train" in logged
    assert "proto_loss_train" in logged

    # get_latent_representation
    latent = model.get_latent_representation()
    assert latent.shape == (adata.n_obs, 5)

    # get_prototypes — shape (n_known_labels, n_latent)
    protos = model.get_prototypes()
    assert protos.shape == (model.n_labels, 5)

    # predict — hard labels
    preds = model.predict()
    assert len(preds) == adata.n_obs
    assert all(p in list(adata.obs["labels"].unique()) for p in preds)

    # predict — subset
    preds_sub = model.predict(indices=[0, 1, 2])
    assert len(preds_sub) == 3

    # predict — soft
    soft = model.predict(soft=True)
    assert isinstance(soft, pd.DataFrame)
    assert soft.shape[0] == adata.n_obs
    assert soft.shape[1] == model.n_labels
    assert np.allclose(soft.values.sum(axis=1), 1.0, atol=1e-5)

    # Inference methods inherited from VAEMixin / RNASeqMixin
    model.get_elbo()
    model.get_reconstruction_error()
    model.get_normalized_expression()

    # Predict on a fresh (but compatible) adata
    adata2 = synthetic_iid()
    preds2 = model.predict(adata2, indices=[0, 1, 2])
    assert len(preds2) == 3


def test_scpoli_all_labeled():
    """All cells are labeled (no cell carries the unlabeled category)."""
    unknown_label = "Unknown"
    adata = synthetic_iid()
    # No cell has label "Unknown" → all cells are labeled
    ScPoli.setup_anndata(
        adata, labels_key="labels", unlabeled_category=unknown_label, batch_key="batch"
    )
    model = ScPoli(adata, n_latent=5)
    assert len(model._unlabeled_indices) == 0
    assert len(model._labeled_indices) == adata.n_obs
    model.train(max_epochs=3, pretrain_epochs=2, train_size=0.9)
    assert model.module._prototypes_initialized.item()


def test_scpoli_all_unlabeled():
    """All cells carry the unlabeled category (no labeled prototypes to train)."""
    adata = synthetic_iid()
    # Overwrite all labels with the unlabeled string
    adata.obs["labels_unlabeled"] = _UNLABELED
    ScPoli.setup_anndata(
        adata,
        labels_key="labels_unlabeled",
        unlabeled_category=_UNLABELED,
        batch_key="batch",
    )
    model = ScPoli(adata, n_latent=5)
    assert model.n_labels == 0
    assert len(model._labeled_indices) == 0
    # Training should not crash even with 0 labeled prototypes
    model.train(max_epochs=2, pretrain_epochs=1, train_size=0.9)


def test_scpoli_multiple_covariates():
    """Continuous and categorical extra covariates are handled correctly."""
    adata = synthetic_iid()
    adata.obs["cont1"] = np.random.normal(size=adata.n_obs)
    adata.obs["cont2"] = np.random.normal(size=adata.n_obs)
    adata.obs["cat1"] = np.random.randint(0, 3, size=adata.n_obs).astype(str)
    adata.obs["cat2"] = np.random.randint(0, 4, size=adata.n_obs).astype(str)

    model, _ = _setup_and_train(
        adata=adata,
        continuous_covariate_keys=["cont1", "cont2"],
        categorical_covariate_keys=["cat1", "cat2"],
    )
    model.get_latent_representation()
    model.get_elbo()
    model.get_normalized_expression()


def test_scpoli_layer():
    """Custom data layer is used instead of adata.X."""
    adata = synthetic_iid()
    adata.layers["counts"] = adata.X.copy()
    model, _ = _setup_and_train(adata=adata, layer="counts")
    model.get_latent_representation()
    model.get_normalized_expression()


def test_scpoli_size_factor():
    """Size-factor key disables the library-size encoder."""
    adata = synthetic_iid()
    adata.obs["size_factor"] = np.random.randint(1, 5, size=adata.n_obs).astype(float)
    model, _ = _setup_and_train(adata=adata, size_factor_key="size_factor")
    model.get_latent_representation()


# ──────────────────────────────────────────────────────────────────────────────
# Prototype-specific tests
# ──────────────────────────────────────────────────────────────────────────────


def test_scpoli_get_prototypes_before_training_raises():
    """get_prototypes() raises if training has not reached prototype phase."""
    adata = synthetic_iid()
    ScPoli.setup_anndata(
        adata, labels_key="labels", unlabeled_category=_UNLABELED, batch_key="batch"
    )
    model = ScPoli(adata, n_latent=5)
    with pytest.raises(RuntimeError, match="Prototypes have not been initialised"):
        model.get_prototypes()


def test_scpoli_get_prototypes_after_training():
    """get_prototypes() returns a DataFrame of shape (n_labels, n_latent) after training."""
    model, adata = _setup_and_train()
    protos = model.get_prototypes()
    assert isinstance(protos, pd.DataFrame)
    assert protos.shape == (model.n_labels, 5)  # n_latent=5 from _setup_and_train default
    # Index should be cell-type label names
    assert list(protos.index) == [model._code_to_label[i] for i in range(model.n_labels)]
    # Columns should be dim_0, dim_1, ...
    assert list(protos.columns) == [f"dim_{i}" for i in range(5)]


def test_scpoli_prototypes_non_trivial():
    """Prototypes should not remain all-zero after training on labeled data."""
    model, _ = _setup_and_train(max_epochs=5, pretrain_epochs=3)
    protos = model.get_prototypes()
    # At least one prototype should be non-zero
    assert not np.allclose(protos, 0.0)


def test_scpoli_pretrain_epochs_zero():
    """pretrain_epochs=0 skips pre-training; prototypes initialized at epoch 0."""
    model, _ = _setup_and_train(max_epochs=3, pretrain_epochs=0)
    assert model.module._prototypes_initialized.item()
    model.get_prototypes()


def test_scpoli_pretrain_epochs_equals_max_epochs():
    """pretrain_epochs == max_epochs: prototype phase never runs."""
    adata = synthetic_iid()
    ScPoli.setup_anndata(
        adata, labels_key="labels", unlabeled_category=_UNLABELED, batch_key="batch"
    )
    model = ScPoli(adata, n_latent=5)
    model.train(max_epochs=3, pretrain_epochs=3, train_size=0.9)
    # Prototype boundary is never crossed → not initialized
    assert not model.module._prototypes_initialized.item()


def test_scpoli_unlabeled_prototypes_kmeans():
    """KMeans-based unlabeled prototype training completes without error."""
    adata = synthetic_iid()
    ScPoli.setup_anndata(
        adata, labels_key="labels", unlabeled_category=_UNLABELED, batch_key="batch"
    )
    model = ScPoli(adata, n_latent=5, unlabeled_weight=0.01)
    model.train(
        max_epochs=4,
        pretrain_epochs=2,
        train_size=0.9,
        unlabeled_prototype_training=True,
        clustering="kmeans",
        n_clusters=4,
    )
    assert model.module._unlabeled_prototypes_initialized.item()
    unlabeled_protos = model.get_unlabeled_prototypes()
    assert unlabeled_protos.shape == (4, 5)


def test_scpoli_unlabeled_prototypes_leiden():
    """Leiden-based unlabeled prototype training completes without error."""
    adata = synthetic_iid()
    ScPoli.setup_anndata(
        adata, labels_key="labels", unlabeled_category=_UNLABELED, batch_key="batch"
    )
    model = ScPoli(adata, n_latent=5, unlabeled_weight=0.01)
    model.train(
        max_epochs=4,
        pretrain_epochs=2,
        train_size=0.9,
        unlabeled_prototype_training=True,
        clustering="leiden",
        clustering_res=1.0,
    )
    assert model.module._unlabeled_prototypes_initialized.item()
    unlabeled_protos = model.get_unlabeled_prototypes()
    assert unlabeled_protos.ndim == 2
    assert unlabeled_protos.shape[1] == 5


def test_scpoli_unlabeled_prototypes_disabled_by_weight():
    """unlabeled_weight=0 (default) means unlabeled prototypes are never seeded."""
    model, _ = _setup_and_train(unlabeled_weight=0.0)
    assert not model.module._unlabeled_prototypes_initialized.item()
    with pytest.raises(RuntimeError, match="Unlabeled prototypes have not been initialised"):
        model.get_unlabeled_prototypes()


# ──────────────────────────────────────────────────────────────────────────────
# Architecture / hyperparameter tests
# ──────────────────────────────────────────────────────────────────────────────


def test_scpoli_batch_embedding_default():
    """Default batch_embedding_dim=10 creates an Embedding layer in the module."""
    adata = synthetic_iid()  # 2 batches
    ScPoli.setup_anndata(
        adata, labels_key="labels", unlabeled_category=_UNLABELED, batch_key="batch"
    )
    model = ScPoli(adata, n_latent=5, batch_embedding_dim=10)
    # The embedding should be registered in the module
    embedding_modules = [
        name for name, _ in model.module.named_modules() if "embedding" in name.lower()
    ]
    assert len(embedding_modules) > 0
    batch_representation = model.get_batch_representation()
    assert isinstance(batch_representation, np.ndarray)
    assert batch_representation.shape == (adata.n_obs, 10)
    model.train(max_epochs=2, pretrain_epochs=1, train_size=0.9)


def test_scpoli_batch_embedding_disabled():
    """batch_embedding_dim=0 disables the Embedding and falls back to one-hot."""
    adata = synthetic_iid()
    ScPoli.setup_anndata(
        adata, labels_key="labels", unlabeled_category=_UNLABELED, batch_key="batch"
    )
    model = ScPoli(adata, n_latent=5, batch_embedding_dim=0)
    model.train(max_epochs=2, pretrain_epochs=1, train_size=0.9)
    model.get_latent_representation()


def test_scpoli_kl_warmup_injected():
    """n_epochs_kl_warmup is injected into plan_kwargs and does not crash."""
    model, _ = _setup_and_train(max_epochs=4, pretrain_epochs=2)
    # TrainingPlan exposes n_epochs_kl_warmup via history indirectly;
    # simplest check: training completed and history exists
    assert "elbo_train" in model.history


def test_scpoli_no_batch_key():
    """Single-batch data (no batch_key) trains without error."""
    adata = synthetic_iid()
    ScPoli.setup_anndata(
        adata,
        labels_key="labels",
        unlabeled_category=_UNLABELED,
        batch_key=None,
    )
    model = ScPoli(adata, n_latent=5, batch_embedding_dim=0)
    model.train(max_epochs=2, pretrain_epochs=1, train_size=0.9)
    model.get_latent_representation()


# ──────────────────────────────────────────────────────────────────────────────
# Query mapping (ArchesMixin.load_query_data)
# ──────────────────────────────────────────────────────────────────────────────


def test_scpoli_online_update(save_path):
    """Reference training + query mapping via load_query_data."""
    n_latent = 5

    # ── Reference ──────────────────────────────────────────────────────────
    adata_ref = synthetic_iid()
    ScPoli.setup_anndata(
        adata_ref, labels_key="labels", unlabeled_category=_UNLABELED, batch_key="batch"
    )
    model_ref = ScPoli(adata_ref, n_latent=n_latent)
    model_ref.train(max_epochs=3, pretrain_epochs=2, train_size=0.9)
    dir_path = os.path.join(save_path, "scpoli_ref/")
    model_ref.save(dir_path, overwrite=True)

    # ── Query — all missing labels ──────────────────────────────────────────
    adata_q = synthetic_iid()
    adata_q.obs["batch"] = adata_q.obs["batch"].cat.rename_categories(["batch_2", "batch_3"])
    adata_q.obs["labels"] = _UNLABELED

    model_q = ScPoli.load_query_data(adata_q, dir_path)
    q_embedding = model_q.module.get_embedding(scvi.REGISTRY_KEYS.BATCH_KEY)
    n_ref_batches = model_ref.module.get_embedding(scvi.REGISTRY_KEYS.BATCH_KEY).num_embeddings
    q_embedding_before = q_embedding.weight.detach().clone()
    protos_before = model_q.module.prototypes_labeled.detach().clone()
    assert q_embedding.weight.requires_grad
    model_q.train(max_epochs=2, pretrain_epochs=0)
    q_embedding_after = q_embedding.weight.detach()
    assert (
        torch.linalg.vector_norm(
            q_embedding_after[n_ref_batches:] - q_embedding_before[n_ref_batches:]
        )
        > 0
    )
    torch.testing.assert_close(model_q.module.prototypes_labeled, protos_before)
    latent_q = model_q.get_latent_representation()
    assert latent_q.shape == (adata_q.n_obs, n_latent)
    preds_q = model_q.predict()
    assert len(preds_q) == adata_q.n_obs

    # ── Query — some labels known ───────────────────────────────────────────
    adata_q2 = synthetic_iid()
    adata_q2.obs["batch"] = adata_q2.obs["batch"].cat.rename_categories(["batch_4", "batch_5"])
    model_q2 = ScPoli.load_query_data(adata_q2, dir_path)
    protos_q2_before = model_q2.module.prototypes_labeled.detach().clone()
    model_q2.train(max_epochs=2, pretrain_epochs=0)
    torch.testing.assert_close(model_q2.module.prototypes_labeled, protos_q2_before)
    model_q2.get_latent_representation()
    model_q2.predict()


# ──────────────────────────────────────────────────────────────────────────────
# Saving / loading
# ──────────────────────────────────────────────────────────────────────────────


def test_scpoli_saving_and_loading(save_path):
    """Saved and reloaded model produces identical predictions."""
    model, adata = _setup_and_train()
    preds_before = model.predict()
    protos_before = model.get_prototypes()

    save_dir = os.path.join(save_path, "scpoli_save/")
    model.save(save_dir, overwrite=True, save_anndata=True)

    # Load with saved adata
    model2 = ScPoli.load(save_dir)
    preds_after = model2.predict()
    protos_after = model2.get_prototypes()

    np.testing.assert_array_equal(preds_before, preds_after)
    np.testing.assert_array_almost_equal(protos_before, protos_after)
    assert model2.is_trained is True

    # Load with user-supplied adata
    model3 = ScPoli.load(save_dir, adata=adata)
    preds3 = model3.predict()
    np.testing.assert_array_equal(preds_before, preds3)

    # Loading with wrong gene set raises
    adata_wrong = synthetic_iid(n_genes=200)
    with pytest.raises(ValueError):
        ScPoli.load(save_dir, adata=adata_wrong)

    # Continue training after loading
    model2.train(max_epochs=2, pretrain_epochs=0)


def test_scpoli_saving_with_prototypes(save_path):
    """Prototype buffers are correctly persisted in the saved state dict."""
    model, _ = _setup_and_train()
    protos = model.get_prototypes()

    save_dir = os.path.join(save_path, "scpoli_proto_save/")
    model.save(save_dir, overwrite=True, save_anndata=True)

    model2 = ScPoli.load(save_dir)
    protos2 = model2.get_prototypes()
    np.testing.assert_array_almost_equal(protos, protos2)

    # The initialized flag is also saved
    assert model2.module._prototypes_initialized.item()


# ──────────────────────────────────────────────────────────────────────────────
# Minified data
# ──────────────────────────────────────────────────────────────────────────────


def _prep_minified_model():
    """Helper: train model, store qzm/qzv, return (model, adata_lib_size)."""
    model, adata = _setup_and_train()
    adata_counts = model.adata_manager.get_from_registry(scvi.REGISTRY_KEYS.X_KEY)
    adata_lib_size = np.squeeze(np.asarray(adata_counts.sum(axis=1)))

    scvi.settings.seed = 1
    qzm, qzv = model.get_latent_representation(give_mean=False, return_dist=True)
    model.adata.obsm["X_latent_qzm"] = qzm
    model.adata.obsm["X_latent_qzv"] = qzv
    return model, adata_lib_size


def test_scpoli_minify_adata():
    """minify_adata() replaces count data with zeros and stores lib size."""
    model, adata_lib_size = _prep_minified_model()
    adata_before = model.adata.copy()

    model.minify_adata()

    assert model.minified_data_type == ADATA_MINIFY_TYPE.LATENT_POSTERIOR
    assert model.adata_manager.registry is model.registry_

    # Original adata object is not modified in place
    assert _is_minified(adata_before) is False

    # Observed lib size is stored correctly
    assert _SCPOLI_OBSERVED_LIB_SIZE in model.adata.obs.columns
    np.testing.assert_array_almost_equal(
        model.adata.obs[_SCPOLI_OBSERVED_LIB_SIZE].values, adata_lib_size
    )

    # var and obs metadata are preserved
    assert model.adata.var_names.equals(adata_before.var_names)


def test_scpoli_minified_latent_representation():
    """Latent representation is consistent before and after minification."""
    model, _ = _prep_minified_model()

    scvi.settings.seed = 1
    latent_orig = model.get_latent_representation()

    model.minify_adata()
    assert model.minified_data_type == ADATA_MINIFY_TYPE.LATENT_POSTERIOR

    scvi.settings.seed = 1
    latent_mini = model.get_latent_representation()

    np.testing.assert_array_equal(latent_mini, latent_orig)


def test_scpoli_minified_predict():
    """predict() works on minified adata and returns consistent labels."""
    model, _ = _prep_minified_model()
    preds_orig = model.predict()

    model.minify_adata()
    preds_mini = model.predict()

    np.testing.assert_array_equal(preds_mini, preds_orig)


def test_scpoli_minified_get_prototypes():
    """get_prototypes() is unaffected by minification (reads module buffer)."""
    model, _ = _prep_minified_model()
    protos_orig = model.get_prototypes()

    model.minify_adata()
    protos_mini = model.get_prototypes()

    np.testing.assert_array_almost_equal(protos_mini, protos_orig)


def test_scpoli_minified_unsupported_methods():
    """Methods requiring count data raise ValueError on minified adata."""
    model, _ = _prep_minified_model()
    model.minify_adata()

    common_msg = "The {} function currently does not support minified data."

    with pytest.raises(ValueError) as exc:
        model.get_elbo()
    assert str(exc.value) == common_msg.format("VAEMixin.get_elbo")

    with pytest.raises(ValueError) as exc:
        model.get_reconstruction_error()
    assert str(exc.value) == common_msg.format("VAEMixin.get_reconstruction_error")

    with pytest.raises(ValueError) as exc:
        model.get_marginal_ll()
    assert str(exc.value) == common_msg.format("VAEMixin.get_marginal_ll")

    with pytest.raises(ValueError) as exc:
        model.get_latent_library_size()
    assert str(exc.value) == common_msg.format("RNASeqMixin.get_latent_library_size")


def test_scpoli_minified_save_then_load(save_path):
    """Minified adata is saved and re-loaded correctly; minified_data_type persists."""
    model, _ = _prep_minified_model()

    scvi.settings.seed = 1
    latent_orig = model.get_latent_representation()

    model.minify_adata()
    save_dir = os.path.join(save_path, "scpoli_mini_save/")
    model.save(save_dir, overwrite=True, save_anndata=True)

    loaded = ScPoli.load(save_dir)
    assert loaded.minified_data_type == ADATA_MINIFY_TYPE.LATENT_POSTERIOR

    scvi.settings.seed = 1
    latent_loaded = loaded.get_latent_representation()
    np.testing.assert_array_equal(latent_loaded, latent_orig)


def test_scpoli_minified_save_then_load_non_minified_adata(save_path):
    """Loading a minified save with a non-minified adata clears minified_data_type."""
    model, adata = _setup_and_train()
    adata_orig = adata.copy()

    qzm, qzv = model.get_latent_representation(give_mean=False, return_dist=True)
    model.adata.obsm["X_latent_qzm"] = qzm
    model.adata.obsm["X_latent_qzv"] = qzv
    model.minify_adata()

    save_dir = os.path.join(save_path, "scpoli_mini_save2/")
    model.save(save_dir, overwrite=True, save_anndata=True)

    loaded = ScPoli.load(save_dir, adata=adata_orig)
    assert loaded.minified_data_type is None
