"""Tests for the low-count model family.

Covers:
- NonZeroSCVI
- ThinnedSCVI
- JointEmbeddingSCVI
- DeterministicThinnedSCVI
"""

from __future__ import annotations

import numpy as np
import pytest
import scipy.sparse as sp

from scvi.data import synthetic_iid
from scvi.external import (
    DeterministicThinnedSCVI,
    JointEmbeddingSCVI,
    NonZeroSCVI,
    ThinnedSCVI,
)

# ---------------------------------------------------------------------------
# NonZeroSCVI
# ---------------------------------------------------------------------------


@pytest.mark.parametrize("n_latent", [5])
def test_nonzero_scvi_basic(n_latent: int):
    """Basic training and inference of NonZeroSCVI."""
    adata = synthetic_iid()
    NonZeroSCVI.setup_anndata(adata, batch_key="batch", labels_key="labels")
    model = NonZeroSCVI(adata, n_latent=n_latent)
    model.train(1, train_size=0.5)
    assert model.is_trained is True
    z = model.get_latent_representation()
    assert z.shape == (adata.shape[0], n_latent)
    model.get_elbo()
    model.get_marginal_ll(n_mc_samples=3)
    model.get_reconstruction_error()
    model.get_normalized_expression()


def test_nonzero_scvi_sparse():
    """NonZeroSCVI with sparse input."""
    adata = synthetic_iid()
    adata.X = sp.csr_matrix(adata.X)
    NonZeroSCVI.setup_anndata(adata)
    model = NonZeroSCVI(adata, n_latent=5)
    model.train(1, train_size=0.5)
    assert model.is_trained is True
    z = model.get_latent_representation()
    assert z.shape == (adata.shape[0], 5)


@pytest.mark.parametrize("normalize_by_nonzero", [True, False])
def test_nonzero_scvi_normalization_modes(normalize_by_nonzero: bool):
    """Both normalization modes for NonZeroSCVI."""
    adata = synthetic_iid()
    NonZeroSCVI.setup_anndata(adata)
    model = NonZeroSCVI(adata, n_latent=5, normalize_by_nonzero=normalize_by_nonzero)
    model.train(1, train_size=0.5)
    assert model.is_trained is True
    z = model.get_latent_representation()
    assert z.shape == (adata.shape[0], 5)


@pytest.mark.parametrize("gene_likelihood", ["zinb", "nb", "poisson"])
def test_nonzero_scvi_gene_likelihoods(gene_likelihood: str):
    """NonZeroSCVI with different gene likelihoods."""
    adata = synthetic_iid()
    NonZeroSCVI.setup_anndata(adata)
    model = NonZeroSCVI(adata, n_latent=5, gene_likelihood=gene_likelihood)
    model.train(1, train_size=0.5)
    assert model.is_trained is True


def test_nonzero_scvi_model_summary():
    """Model summary string includes normalize_by_nonzero."""
    adata = synthetic_iid()
    NonZeroSCVI.setup_anndata(adata)

    model_true = NonZeroSCVI(adata, normalize_by_nonzero=True)
    assert "normalize_by_nonzero: True" in model_true._model_summary_string

    model_false = NonZeroSCVI(adata, normalize_by_nonzero=False)
    assert "normalize_by_nonzero: False" in model_false._model_summary_string


def test_nonzero_scvi_masking_effect():
    """NonZeroSCVI trains correctly when some genes are all-zero."""
    adata = synthetic_iid()
    adata.X[:, :10] = 0
    NonZeroSCVI.setup_anndata(adata)
    model = NonZeroSCVI(adata, n_latent=5)
    model.train(1, train_size=0.5)
    assert model.is_trained is True


# ---------------------------------------------------------------------------
# ThinnedSCVI
# ---------------------------------------------------------------------------


@pytest.mark.parametrize("n_latent", [5])
def test_thinned_scvi_basic(n_latent: int):
    """Basic training and inference of ThinnedSCVI."""
    adata = synthetic_iid()
    ThinnedSCVI.setup_anndata(adata, batch_key="batch", labels_key="labels")
    model = ThinnedSCVI(adata, n_latent=n_latent)
    model.train(1, train_size=0.5)
    assert model.is_trained is True
    z = model.get_latent_representation()
    assert z.shape == (adata.shape[0], n_latent)
    model.get_elbo()
    model.get_reconstruction_error()
    model.get_normalized_expression()


def test_thinned_scvi_sparse():
    """ThinnedSCVI with sparse input."""
    adata = synthetic_iid()
    adata.X = sp.csr_matrix(adata.X)
    ThinnedSCVI.setup_anndata(adata)
    model = ThinnedSCVI(adata, n_latent=5)
    model.train(1, train_size=0.5)
    assert model.is_trained is True
    z = model.get_latent_representation()
    assert z.shape == (adata.shape[0], 5)


@pytest.mark.parametrize("min_library_size", [10.0, 50.0, 100.0])
def test_thinned_scvi_min_library_size(min_library_size: float):
    """ThinnedSCVI with different min_library_size values."""
    adata = synthetic_iid()
    ThinnedSCVI.setup_anndata(adata)
    model = ThinnedSCVI(adata, n_latent=5, min_library_size=min_library_size)
    model.train(1, train_size=0.5)
    assert model.is_trained is True


def test_thinned_scvi_model_summary():
    """Model summary string includes min_library_size."""
    adata = synthetic_iid()
    ThinnedSCVI.setup_anndata(adata)
    model = ThinnedSCVI(adata, min_library_size=42.0)
    assert "min_library_size: 42.0" in model._model_summary_string


def test_thinned_scvi_inference_uses_full_data():
    """At inference time, ThinnedSCVI returns identical z on two calls."""
    adata = synthetic_iid()
    ThinnedSCVI.setup_anndata(adata)
    model = ThinnedSCVI(adata, n_latent=5)
    model.train(1, train_size=0.5)

    z1 = model.get_latent_representation()
    z2 = model.get_latent_representation()
    np.testing.assert_array_equal(z1, z2)


# ---------------------------------------------------------------------------
# JointEmbeddingSCVI
# ---------------------------------------------------------------------------


@pytest.mark.parametrize("n_latent", [5])
def test_joint_embedding_scvi_basic(n_latent: int):
    """Basic training and inference of JointEmbeddingSCVI."""
    adata = synthetic_iid()
    JointEmbeddingSCVI.setup_anndata(adata, batch_key="batch", labels_key="labels")
    model = JointEmbeddingSCVI(adata, n_latent=n_latent)
    model.train(1, train_size=0.5)
    assert model.is_trained is True
    z = model.get_latent_representation()
    assert z.shape == (adata.shape[0], n_latent)
    model.get_elbo()
    model.get_reconstruction_error()
    model.get_normalized_expression()


def test_joint_embedding_scvi_sparse():
    """JointEmbeddingSCVI with sparse input."""
    adata = synthetic_iid()
    adata.X = sp.csr_matrix(adata.X)
    JointEmbeddingSCVI.setup_anndata(adata)
    model = JointEmbeddingSCVI(adata, n_latent=5)
    model.train(1, train_size=0.5)
    assert model.is_trained is True
    z = model.get_latent_representation()
    assert z.shape == (adata.shape[0], 5)


@pytest.mark.parametrize(
    ("joint_embedding_weight", "reconstruction_weight"),
    [(1.0, 1.0), (0.5, 1.0), (1.0, 0.0), (0.0, 1.0)],
)
def test_joint_embedding_scvi_loss_weights(
    joint_embedding_weight: float, reconstruction_weight: float
):
    """JointEmbeddingSCVI with various loss weight combinations."""
    adata = synthetic_iid()
    JointEmbeddingSCVI.setup_anndata(adata)
    model = JointEmbeddingSCVI(
        adata,
        n_latent=5,
        joint_embedding_weight=joint_embedding_weight,
        reconstruction_weight=reconstruction_weight,
    )
    model.train(1, train_size=0.5)
    assert model.is_trained is True


def test_joint_embedding_scvi_disabled():
    """JointEmbeddingSCVI with use_joint_embedding=False behaves like SCVI."""
    adata = synthetic_iid()
    JointEmbeddingSCVI.setup_anndata(adata)
    model = JointEmbeddingSCVI(adata, n_latent=5, use_joint_embedding=False)
    model.train(1, train_size=0.5)
    assert model.is_trained is True
    z = model.get_latent_representation()
    assert z.shape == (adata.shape[0], 5)


def test_joint_embedding_scvi_variance_weight():
    """JointEmbeddingSCVI with variance regularization enabled."""
    adata = synthetic_iid()
    JointEmbeddingSCVI.setup_anndata(adata)
    model = JointEmbeddingSCVI(adata, n_latent=5, variance_weight=1.0)
    model.train(1, train_size=0.5)
    assert model.is_trained is True


def test_joint_embedding_scvi_model_summary():
    """Model summary string contains key hyperparameters."""
    adata = synthetic_iid()
    JointEmbeddingSCVI.setup_anndata(adata)
    model = JointEmbeddingSCVI(
        adata,
        joint_embedding_weight=2.0,
        lambda_off_diag=0.05,
        min_library_size=20.0,
        reconstruction_weight=0.5,
        variance_weight=1.0,
    )
    s = model._model_summary_string
    assert "joint_embedding_weight: 2.0" in s
    assert "lambda_off_diag: 0.05" in s
    assert "min_library_size: 20.0" in s
    assert "reconstruction_weight: 0.5" in s
    assert "variance_weight: 1.0" in s


@pytest.mark.parametrize("gene_likelihood", ["zinb", "nb", "poisson"])
def test_joint_embedding_scvi_gene_likelihoods(gene_likelihood: str):
    """JointEmbeddingSCVI with different gene likelihoods."""
    adata = synthetic_iid()
    JointEmbeddingSCVI.setup_anndata(adata)
    model = JointEmbeddingSCVI(adata, n_latent=5, gene_likelihood=gene_likelihood)
    model.train(1, train_size=0.5)
    assert model.is_trained is True


# ---------------------------------------------------------------------------
# DeterministicThinnedSCVI
# ---------------------------------------------------------------------------


@pytest.mark.parametrize("n_latent", [5])
def test_deterministic_thinned_scvi_basic(n_latent: int):
    """Basic training and inference of DeterministicThinnedSCVI."""
    adata = synthetic_iid()
    DeterministicThinnedSCVI.setup_anndata(adata, batch_key="batch", labels_key="labels")
    model = DeterministicThinnedSCVI(adata, n_latent=n_latent)
    model.train(1, train_size=0.5)
    assert model.is_trained is True
    z = model.get_latent_representation()
    assert z.shape == (adata.shape[0], n_latent)
    model.get_elbo()
    model.get_reconstruction_error()
    model.get_normalized_expression()


def test_deterministic_thinned_scvi_sparse():
    """DeterministicThinnedSCVI with sparse input."""
    adata = synthetic_iid()
    adata.X = sp.csr_matrix(adata.X)
    DeterministicThinnedSCVI.setup_anndata(adata)
    model = DeterministicThinnedSCVI(adata, n_latent=5)
    model.train(1, train_size=0.5)
    assert model.is_trained is True
    z = model.get_latent_representation()
    assert z.shape == (adata.shape[0], 5)


def test_deterministic_thinned_scvi_model_summary():
    """Model summary string includes min_library_size."""
    adata = synthetic_iid()
    DeterministicThinnedSCVI.setup_anndata(adata)
    model = DeterministicThinnedSCVI(adata, min_library_size=42.0)
    assert "min_library_size: 42.0" in model._model_summary_string


def test_deterministic_thinned_scvi_deterministic_latent():
    """At inference time, DeterministicThinnedSCVI returns identical z on two calls."""
    adata = synthetic_iid()
    DeterministicThinnedSCVI.setup_anndata(adata)
    model = DeterministicThinnedSCVI(adata, n_latent=5)
    model.train(1, train_size=0.5)

    z1 = model.get_latent_representation()
    z2 = model.get_latent_representation()
    np.testing.assert_array_equal(z1, z2)


@pytest.mark.parametrize("min_library_size", [10.0, 50.0])
def test_deterministic_thinned_scvi_min_library_size(min_library_size: float):
    """DeterministicThinnedSCVI with different min_library_size values."""
    adata = synthetic_iid()
    DeterministicThinnedSCVI.setup_anndata(adata)
    model = DeterministicThinnedSCVI(adata, n_latent=5, min_library_size=min_library_size)
    model.train(1, train_size=0.5)
    assert model.is_trained is True
