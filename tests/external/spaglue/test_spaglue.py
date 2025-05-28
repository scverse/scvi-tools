import numpy as np
import pandas as pd
import pytest
import torch
from anndata import AnnData

from scvi.external.spaglue._model import SPAGLUE


# Fixtures for creating test data
@pytest.fixture
def adata_seq():
    """Create a mock AnnData object for sequencing data."""
    n_obs = 100  # Number of cells
    n_vars = 50  # Number of genes
    X = np.random.poisson(1.0, size=(n_obs, n_vars))
    obs = {"batch": pd.Categorical(np.random.choice(["batch1", "batch2"], size=n_obs))}
    return AnnData(X=X, obs=obs)


@pytest.fixture
def adata_spatial():
    """Create a mock AnnData object for spatial data."""
    n_obs = 80  # Number of spots
    n_vars = 40  # Number of genes
    X = np.random.poisson(1.0, size=(n_obs, n_vars))
    obs = {"batch": pd.Categorical(np.random.choice(["batch1", "batch2"], size=n_obs))}
    return AnnData(X=X, obs=obs)


# Test initialization
def test_spaglue_initialization(adata_seq, adata_spatial):
    """Test the initialization of the SPAGLUE model."""
    SPAGLUE.setup_anndata(adata_seq, batch_key="batch")
    SPAGLUE.setup_anndata(adata_spatial, batch_key="batch")

    model = SPAGLUE(adata_seq, adata_spatial)
    assert model.module is not None
    assert model.adatas[0] is adata_seq
    assert model.adatas[1] is adata_spatial


# Test training
def test_spaglue_training(adata_seq, adata_spatial):
    """Test the training process of the SPAGLUE model."""
    SPAGLUE.setup_anndata(adata_seq, batch_key="batch")
    SPAGLUE.setup_anndata(adata_spatial, batch_key="batch")

    model = SPAGLUE(adata_seq, adata_spatial)
    model.train(max_epochs=1, batch_size=16)

    assert model.is_trained_ is True


# Test latent representation
def test_spaglue_latent_representation(adata_seq, adata_spatial):
    """Test the computation of latent representations."""
    SPAGLUE.setup_anndata(adata_seq, batch_key="batch")
    SPAGLUE.setup_anndata(adata_spatial, batch_key="batch")

    model = SPAGLUE(adata_seq, adata_spatial)
    model.train(max_epochs=1, batch_size=16)

    latent = model.get_latent_representation(adata_seq=adata_seq, adata_spatial=adata_spatial)
    assert "seq" in latent
    assert "spatial" in latent
    assert latent["seq"].shape[0] == adata_seq.shape[0]
    assert latent["spatial"].shape[0] == adata_spatial.shape[0]


# Test generative model
def test_spaglue_generative_model(adata_seq, adata_spatial):
    """Test the generative model outputs."""
    SPAGLUE.setup_anndata(adata_seq, batch_key="batch")
    SPAGLUE.setup_anndata(adata_spatial, batch_key="batch")

    model = SPAGLUE(adata_seq, adata_spatial)
    model.train(max_epochs=1, batch_size=16)

    # Test generative outputs
    tensors = {
        "X": torch.tensor(adata_seq.X, dtype=torch.float32),
        "batch": torch.tensor(adata_seq.obs["batch"].cat.codes.values, dtype=torch.long),
    }
    inference_outputs = model.module.inference(tensors["X"], mode=0)
    generative_outputs = model.module.generative(
        inference_outputs["z"],
        inference_outputs["library"],
        batch_index=tensors["batch"],
        mode=0,
    )

    assert "px" in generative_outputs
    assert "pz" in generative_outputs


# Test loss computation
def test_spaglue_loss(adata_seq, adata_spatial):
    """Test the loss computation."""
    SPAGLUE.setup_anndata(adata_seq, batch_key="batch")
    SPAGLUE.setup_anndata(adata_spatial, batch_key="batch")

    model = SPAGLUE(adata_seq, adata_spatial)
    model.train(max_epochs=1, batch_size=16)

    device = next(model.module.parameters()).device

    tensors = {
        "X": torch.tensor(adata_seq.X, dtype=torch.float32, device=device),
        "batch": torch.tensor(
            adata_seq.obs["batch"].cat.codes.values, dtype=torch.long, device=device
        ),
    }
    inference_outputs = model.module.inference(tensors["X"], mode=0)
    generative_outputs = model.module.generative(
        inference_outputs["z"],
        inference_outputs["library"],
        batch_index=tensors["batch"],
        mode=0,
    )
    loss = model.module.loss(tensors, inference_outputs, generative_outputs)

    assert loss.loss is not None
    assert loss.reconstruction_loss is not None
    assert loss.kl_local is not None
