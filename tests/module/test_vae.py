import pytest
import torch

from scvi import REGISTRY_KEYS
from scvi.module import VAE, gaussian_kernel


@pytest.mark.parametrize("n_samples", [1, 2, 3])
@pytest.mark.parametrize("gene_likelihood", ["zinb", "nb", "poisson"])
def test_sample(
    n_samples: int,
    gene_likelihood: str,
    n_input: int = 100,
    batch_size: int = 10,
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


def test_gaussian_kernel():
    # Define sample data
    x = torch.tensor([[1.0, 2.0], [2.0, 3.0], [3.0, 4.0]])
    y = torch.tensor([[2.0, 1.0], [3.0, 2.0]])

    # Compute the Gaussian kernel matrix
    kernel_matrix = gaussian_kernel(x, y)
    print(f"{kernel_matrix = }")

    # Define the expected kernel matrix manually based on the sample data
    expected_kernel_matrix = torch.tensor([[0.1353, 0.0183], [0.0183, 0.1353], [0.0001, 0.0183]])
    print(f"{expected_kernel_matrix = }")

    # Check if the computed kernel matrix matches the expected kernel matrix
    assert torch.allclose(
        kernel_matrix, expected_kernel_matrix, atol=1e-4
    ), "Kernel matrix computation incorrect"
