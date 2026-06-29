import numpy as np
import pytest
import scipy.sparse as sp

from scvi.data import synthetic_iid
from scvi.external import JointEmbeddingSCVI

CCO_METRICS = {"cco_loss", "cco_invariance", "cco_redundancy", "variance_loss"}

N_LATENT = 5


def _get_adata(sparse: bool = False):
    adata = synthetic_iid()
    if sparse:
        adata.X = sp.csr_matrix(adata.X)
    JointEmbeddingSCVI.setup_anndata(adata, batch_key="batch", labels_key="labels")
    return adata


def test_joint_embedding_scvi_train_and_latent():
    adata = _get_adata()
    model = JointEmbeddingSCVI(adata, n_latent=N_LATENT)
    model.train(1, train_size=0.5)
    assert model.is_trained is True
    latent = model.get_latent_representation()
    assert latent.shape == (adata.n_obs, N_LATENT)
    assert np.all(np.isfinite(latent))


def test_joint_embedding_scvi_inference_methods():
    adata = _get_adata()
    model = JointEmbeddingSCVI(adata, n_latent=N_LATENT)
    model.train(1, train_size=0.5)

    assert model.get_elbo().ndim == 0
    assert model.get_marginal_ll(n_mc_samples=3).ndim == 0
    assert model.get_reconstruction_error()["reconstruction_loss"].ndim == 0
    assert model.get_normalized_expression().shape == (adata.n_obs, adata.n_vars)


def test_joint_embedding_scvi_sparse():
    adata = _get_adata(sparse=True)
    model = JointEmbeddingSCVI(adata, n_latent=N_LATENT)
    model.train(1, train_size=0.5)
    assert model.is_trained is True
    assert model.get_latent_representation().shape == (adata.n_obs, N_LATENT)


@pytest.mark.parametrize("gene_likelihood", ["zinb", "nb", "poisson"])
def test_joint_embedding_scvi_gene_likelihood(gene_likelihood):
    adata = _get_adata()
    model = JointEmbeddingSCVI(adata, n_latent=N_LATENT, gene_likelihood=gene_likelihood)
    model.train(1, train_size=0.5)
    assert model.is_trained is True
    assert model.get_latent_representation().shape == (adata.n_obs, N_LATENT)


def test_joint_embedding_scvi_fallback():
    adata = _get_adata()
    model = JointEmbeddingSCVI(adata, n_latent=N_LATENT, use_joint_embedding=False)
    model.train(1, train_size=0.5)
    assert model.is_trained is True
    assert model.get_latent_representation().shape == (adata.n_obs, N_LATENT)


def _one_batch(model, adata):
    """Return one correctly-formatted minibatch of registered tensors on the model's device."""
    dl = model._make_data_loader(adata=adata, indices=np.arange(adata.n_obs))
    tensors = next(iter(dl))
    device = next(model.module.parameters()).device
    return {k: v.to(device) for k, v in tensors.items()}


def test_joint_embedding_scvi_cco_branch_active_in_training():
    """The CCO loss branch only runs in training mode and only when enabled."""
    adata = _get_adata()
    model = JointEmbeddingSCVI(adata, n_latent=N_LATENT)
    model.train(1, train_size=0.5)
    tensors = _one_batch(model, adata)
    module = model.module

    # Training mode + enabled: CCO metrics present and finite
    module.train()
    _, _, losses = module(tensors)
    assert CCO_METRICS.issubset(losses.extra_metrics)
    assert all(np.isfinite(losses.extra_metrics[k].item()) for k in CCO_METRICS)

    # Eval mode: CCO branch is skipped (ELBO only)
    module.eval()
    _, _, losses_eval = module(tensors)
    assert not CCO_METRICS.intersection(losses_eval.extra_metrics)


def test_joint_embedding_scvi_cco_branch_disabled():
    """With use_joint_embedding=False, the CCO branch never runs, even in training."""
    adata = _get_adata()
    model = JointEmbeddingSCVI(adata, n_latent=N_LATENT, use_joint_embedding=False)
    model.train(1, train_size=0.5)
    tensors = _one_batch(model, adata)
    model.module.train()
    _, _, losses = model.module(tensors)
    assert not CCO_METRICS.intersection(losses.extra_metrics)


def test_joint_embedding_scvi_hyperparameters_in_summary():
    adata = _get_adata()
    model = JointEmbeddingSCVI(
        adata,
        n_latent=N_LATENT,
        joint_embedding_weight=2.0,
        lambda_off_diag=0.05,
        min_library_size=20.0,
        reconstruction_weight=0.5,
        variance_weight=1.0,
    )
    summary = model._model_summary_string
    assert "joint_embedding_weight: 2.0" in summary
    assert "lambda_off_diag: 0.05" in summary
    assert "min_library_size: 20.0" in summary
    assert "reconstruction_weight: 0.5" in summary
    assert "variance_weight: 1.0" in summary


def test_joint_embedding_scvi_save_load(tmp_path):
    adata = _get_adata()
    model = JointEmbeddingSCVI(
        adata,
        n_latent=N_LATENT,
        joint_embedding_weight=2.0,
        lambda_off_diag=0.05,
        min_library_size=20.0,
        reconstruction_weight=0.5,
        variance_weight=1.0,
    )
    model.train(1, train_size=0.5)
    latent_before = model.get_latent_representation()

    save_dir = str(tmp_path)
    model.save(save_dir, overwrite=True, save_anndata=True)

    loaded = JointEmbeddingSCVI.load(save_dir)
    latent_after = loaded.get_latent_representation()

    assert np.allclose(latent_before, latent_after)
    # New init params must round-trip through save/load
    assert loaded.module.joint_embedding_weight == 2.0
    assert loaded.module.lambda_off_diag == 0.05
    assert loaded.module.min_library_size == 20.0
    assert loaded.module.reconstruction_weight == 0.5
    assert loaded.module.variance_weight == 1.0
