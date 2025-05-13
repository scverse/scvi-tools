import pytest
import torch

from scvi.data import synthetic_iid
from scvi.external.spaglue import SPAGLUEVAE


@pytest.fixture(scope="session")
def mock_adata():
    """Mock AnnData for SPAGLUE testing."""
    adata = synthetic_iid()
    adata.obs["batch"] = ["batch_0"] * (adata.n_obs // 2) + ["batch_1"] * (adata.n_obs // 2)
    return adata


@pytest.fixture(scope="session")
def spaglue_model(mock_adata):
    """Initialize SPAGLUE model."""
    n_inputs = [mock_adata.n_vars, mock_adata.n_vars]
    gene_likelihoods = ["nb", "nb"]
    model = SPAGLUEVAE(n_inputs=n_inputs, gene_likelihoods=gene_likelihoods, dropout_rate=0.1)
    return model


def test_spaglue_latent_representation(mock_adata, spaglue_model):
    """Test SPAGLUE latent representation extraction."""
    tensors = {"X": torch.tensor(mock_adata.X, dtype=torch.float32)}
    inference_outputs = spaglue_model.inference(tensors["X"], mode=0)
    z = inference_outputs["z"]

    assert z.shape[0] == mock_adata.n_obs, (
        "Latent representation should match number of observations."
    )


def test_spaglue_reconstruction_loss(mock_adata, spaglue_model):
    """Test SPAGLUE reconstruction loss computation."""
    tensors = {"X": torch.tensor(mock_adata.X, dtype=torch.float32)}
    inference_outputs = spaglue_model.inference(tensors["X"], mode=0)
    generative_outputs = spaglue_model.generative(
        z=inference_outputs["z"],
        library=inference_outputs["library"],
        mode=0,
    )
    reconstruction_loss = spaglue_model.reconstruction_loss(
        tensors["X"],
        generative_outputs["mu"],
        generative_outputs["log_theta"],
        mode=0,
    )

    assert reconstruction_loss.shape[0] == mock_adata.n_obs, (
        "Reconstruction loss should match number of observations."
    )
    assert torch.all(reconstruction_loss >= 0), "Reconstruction loss should be non-negative."


def test_spaglue_kl_divergence(mock_adata, spaglue_model):
    """Test SPAGLUE KL divergence computation."""
    tensors = {"X": torch.tensor(mock_adata.X, dtype=torch.float32)}

    inference_outputs = spaglue_model.inference(tensors["X"], mode=0)
    generative_outputs = spaglue_model.generative(
        z=inference_outputs["z"],
        library=inference_outputs["library"],
        mode=0,
    )

    loss_output = spaglue_model.loss(
        tensors,
        inference_outputs,
        generative_outputs,
        mode=0,
    )
    assert loss_output.kl_local["kl_local"].shape[0] == mock_adata.n_obs, (
        "KL divergence should match number of observations."
    )
    assert torch.all(loss_output.kl_local["kl_local"] >= 0), (
        "KL divergence should be non-negative."
    )


def test_spaglue_latent_representation_consistency(mock_adata, spaglue_model):
    """Test that SPAGLUE latent representations are consistent."""
    tensors = {"X": torch.tensor(mock_adata.X, dtype=torch.float32)}
    inference_outputs_1 = spaglue_model.inference(tensors["X"], mode=0)
    inference_outputs_2 = spaglue_model.inference(tensors["X"], mode=0)

    print(inference_outputs_1["qm"].shape)
    print(inference_outputs_1["z"].shape)

    # check that means/var are different due to randomization
    assert not torch.allclose(inference_outputs_1["qm"], inference_outputs_2["qm"], atol=1e-3), (
        "Mean (qm) should differ when dropout rate is applied."
    )
    assert not torch.allclose(inference_outputs_1["qv"], inference_outputs_2["qv"], atol=1e-6), (
        "Variance (qv) should differ when dropout rate is applied."
    )
