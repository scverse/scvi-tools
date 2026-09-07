import pytest
import torch

import scvi.module._vae as vae_module
from scvi import REGISTRY_KEYS
from scvi.module import VAE


@pytest.mark.parametrize("n_samples", [1, 2, 3])
@pytest.mark.parametrize("gene_likelihood", ["zinb", "nb", "poisson"])
@pytest.mark.parametrize("n_input", [100])
@pytest.mark.parametrize("batch_size", [10])
def test_sample(
    n_samples: int,
    gene_likelihood: str,
    n_input: int,
    batch_size: int,
):
    vae = VAE(n_input=n_input, gene_likelihood=gene_likelihood)
    x = torch.randint(0, 100, (batch_size, n_input), dtype=torch.float32)
    batch = torch.zeros(batch_size, dtype=torch.long)
    labels = torch.zeros(batch_size, dtype=torch.long)
    tensors = {
        REGISTRY_KEYS.X_KEY: x,
        REGISTRY_KEYS.BATCH_KEY: batch,
        REGISTRY_KEYS.LABELS_KEY: labels,
    }

    x_hat = vae.sample(tensors, n_samples=n_samples)
    assert x_hat.dtype == torch.float32
    if n_samples > 1:
        assert x_hat.shape == (batch_size, n_input, n_samples)
    else:
        assert x_hat.shape == (batch_size, n_input)


def test_sample_poisson_cpu_detour_forced(monkeypatch):
    """Forces the CPU detour the way a torch without the 'aten::poisson' MPS kernel would."""
    monkeypatch.setattr(vae_module, "_needs_cpu_detour", lambda on_mps, op: True)

    n_input, batch_size = 100, 10
    vae = VAE(n_input=n_input, gene_likelihood="poisson")
    x = torch.randint(0, 100, (batch_size, n_input), dtype=torch.float32)
    batch = torch.zeros(batch_size, dtype=torch.long)
    labels = torch.zeros(batch_size, dtype=torch.long)
    tensors = {
        REGISTRY_KEYS.X_KEY: x,
        REGISTRY_KEYS.BATCH_KEY: batch,
        REGISTRY_KEYS.LABELS_KEY: labels,
    }

    x_hat = vae.sample(tensors, n_samples=1)
    assert x_hat.shape == (batch_size, n_input)
    assert (x_hat >= 0).all()


@pytest.mark.skipif(not torch.backends.mps.is_available(), reason="requires an MPS device")
def test_sample_poisson_stays_on_mps():
    n_input, batch_size = 100, 10
    vae = VAE(n_input=n_input, gene_likelihood="poisson").to("mps")
    x = torch.randint(0, 100, (batch_size, n_input), dtype=torch.float32, device="mps")
    batch = torch.zeros(batch_size, dtype=torch.long, device="mps")
    labels = torch.zeros(batch_size, dtype=torch.long, device="mps")
    tensors = {
        REGISTRY_KEYS.X_KEY: x,
        REGISTRY_KEYS.BATCH_KEY: batch,
        REGISTRY_KEYS.LABELS_KEY: labels,
    }

    x_hat = vae.sample(tensors, n_samples=1)
    assert x_hat.shape == (batch_size, n_input)
