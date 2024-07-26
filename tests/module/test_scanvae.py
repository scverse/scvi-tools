import pytest
import torch

from scvi import REGISTRY_KEYS
from scvi.module import SCANVAE


@pytest.mark.parametrize("n_labels", [1, 2, 10])
def test_scanvae_loss_n_labels(n_labels: int):
    vae = SCANVAE(n_input=100, n_labels=n_labels)

    x = torch.randint(0, 100, (128, 100), dtype=torch.float32)
    batch = torch.randint(0, 5, (128, 1), dtype=torch.long)
    labels = torch.randint(0, 15, (128, 1), dtype=torch.long)
    tensors = {
        REGISTRY_KEYS.X_KEY: x,
        REGISTRY_KEYS.BATCH_KEY: batch,
        REGISTRY_KEYS.LABELS_KEY: labels,
    }

    _ = vae.forward(tensors, inference_kwargs={"n_samples": 1})


@pytest.mark.parametrize("n_samples", [1, 2, 3])
def test_scanvae_loss_n_samples(n_samples: int):
    vae = SCANVAE(n_input=100, n_labels=1)
    vae2 = SCANVAE(n_input=100, n_labels=2)
    vae3 = SCANVAE(n_input=100, n_labels=8)

    x = torch.randint(0, 100, (128, 100), dtype=torch.float32)
    batch = torch.randint(0, 5, (128, 1), dtype=torch.long)
    labels = torch.randint(0, 15, (128, 1), dtype=torch.long)
    tensors = {
        REGISTRY_KEYS.X_KEY: x,
        REGISTRY_KEYS.BATCH_KEY: batch,
        REGISTRY_KEYS.LABELS_KEY: labels,
    }

    _ = vae.forward(tensors, inference_kwargs={"n_samples": n_samples})
    _ = vae2.forward(tensors, inference_kwargs={"n_samples": n_samples})
    _ = vae3.forward(tensors, inference_kwargs={"n_samples": n_samples})
