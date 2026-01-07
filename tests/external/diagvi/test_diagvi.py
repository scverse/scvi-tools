import numpy as np
import pandas as pd
import pytest
import torch
from anndata import AnnData

from scvi.external.diagvi._model import DIAGVI


@pytest.fixture
def adata_seq():
    n_obs = 100
    n_vars = 50
    X = np.random.poisson(1.0, size=(n_obs, n_vars))
    obs = {"batch": pd.Categorical(np.random.choice(["batch1", "batch2"], size=n_obs))}
    var = pd.DataFrame(index=[f"gene{i}" for i in range(n_vars)])
    return AnnData(X=X, obs=obs, var=var)


@pytest.fixture
def adata_spatial():
    n_obs = 80
    n_vars = 50  # Match n_vars and gene names for compatibility
    X = np.random.poisson(1.0, size=(n_obs, n_vars))
    obs = {"batch": pd.Categorical(np.random.choice(["batch1", "batch2"], size=n_obs))}
    var = pd.DataFrame(index=[f"gene{i}" for i in range(n_vars)])
    return AnnData(X=X, obs=obs, var=var)


def make_model(adata_seq, adata_spatial):
    DIAGVI.setup_anndata(adata_seq, batch_key="batch", likelihood="nb")
    DIAGVI.setup_anndata(adata_spatial, batch_key="batch", likelihood="nb")
    return DIAGVI({"diss": adata_seq, "spatial": adata_spatial})


def test_diagvi_initialization(adata_seq, adata_spatial):
    model = make_model(adata_seq, adata_spatial)
    assert model.module is not None
    assert "diss" in model.adatas
    assert "spatial" in model.adatas


def test_diagvi_training(adata_seq, adata_spatial):
    model = make_model(adata_seq, adata_spatial)
    model.train(max_epochs=1, batch_size=16)
    assert model.is_trained_ is True


def test_diagvi_latent_representation(adata_seq, adata_spatial):
    model = make_model(adata_seq, adata_spatial)
    model.train(max_epochs=1, batch_size=16)
    latents = model.get_latent_representation()
    assert latents["diss"].shape[0] == adata_seq.shape[0]
    assert latents["spatial"].shape[0] == adata_spatial.shape[0]


def test_diagvi_generative_model(adata_seq, adata_spatial):
    model = make_model(adata_seq, adata_spatial)
    model.train(max_epochs=1, batch_size=16)
    tensors = {
        "X": torch.tensor(adata_seq.X, dtype=torch.float32),
        "batch": torch.tensor(adata_seq.obs["batch"].cat.codes.values, dtype=torch.long),
    }
    inference_outputs = model.module.inference(tensors["X"], mode="diss")
    generative_outputs = model.module.generative(
        inference_outputs["z"],
        inference_outputs["library"],
        batch_index=tensors["batch"],
        v=inference_outputs["v"],
        mode="diss",
    )
    assert "px" in generative_outputs or "px_rate" in generative_outputs
    assert "pz" in generative_outputs


def test_diagvi_loss(adata_seq, adata_spatial):
    model = make_model(adata_seq, adata_spatial)
    model.train(max_epochs=1, batch_size=16)
    device = next(model.module.parameters()).device
    tensors = {
        "X": torch.tensor(adata_seq.X, dtype=torch.float32, device=device),
        "batch": torch.tensor(
            adata_seq.obs["batch"].cat.codes.values, dtype=torch.long, device=device
        ),
    }
    inference_outputs = model.module.inference(tensors["X"], mode="diss")
    generative_outputs = model.module.generative(
        inference_outputs["z"],
        inference_outputs["library"],
        batch_index=tensors["batch"],
        v=inference_outputs["v"],
        mode="diss",
    )
    loss = model.module.loss(tensors, inference_outputs, generative_outputs, mode="diss")
    assert loss.loss is not None
    assert loss.reconstruction_loss is not None
    assert loss.kl_local is not None


def test_diagvi_get_imputed_values(adata_seq, adata_spatial):
    model = make_model(adata_seq, adata_spatial)
    model.train(max_epochs=1, batch_size=16)
    # Impute spatial from diss
    imputed, _ = model.get_imputed_values(source_name="spatial")
    assert isinstance(imputed, np.ndarray)
    assert imputed.shape[0] == adata_spatial.shape[0]
    assert imputed.shape[1] == adata_seq.shape[1]
    # Impute diss from spatial
    imputed2, _ = model.get_imputed_values(source_name="diss")
    assert isinstance(imputed2, np.ndarray)
    assert imputed2.shape[0] == adata_seq.shape[0]
    assert imputed2.shape[1] == adata_spatial.shape[1]


def test_diagvi_save_and_load(adata_seq, adata_spatial, tmp_path):
    model = make_model(adata_seq, adata_spatial)
    model.train(max_epochs=1, batch_size=16)
    save_path = tmp_path / "diagvi_test_model.pt"
    model.save(save_path, save_anndata=True)

    # Load the model
    loaded_model = DIAGVI.load(save_path)

    # compare state dicts
    model_state = model.module.state_dict()
    loaded_state = loaded_model.module.state_dict()
    for key in model_state:
        assert key in loaded_state
        assert torch.allclose(model_state[key], loaded_state[key], atol=1e-6), f"Mismatch in {key}"

    # compare output equivalence
    latents = loaded_model.get_latent_representation()
    assert latents["diss"].shape[0] == adata_seq.shape[0]
    assert latents["spatial"].shape[0] == adata_spatial.shape[0]


def test_feature_confidence(adata_seq, adata_spatial):
    methods = ["min", "mean", "max", "median"]
    for m in methods:
        model = make_model(adata_seq, adata_spatial)
        # Create a dummy feature embedding (e.g., from PCA or random)
        feature_embedding = np.random.rand(adata_seq.shape[1], 8)
        score = model.compute_per_feature_confidence(feature_embedding, conf_method=m)
        assert isinstance(score, np.ndarray) or isinstance(score, list)
        assert len(score) == feature_embedding.shape[0]
        assert np.all(np.isfinite(score))
