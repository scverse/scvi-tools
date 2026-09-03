import pytest
import torch

from scvi import REGISTRY_KEYS
from scvi.module import VAE
from scvi.module._constants import MODULE_KEYS


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


@pytest.mark.parametrize("pre_encoder_covariate", ["continuous", "embedding_batch"])
def test_regular_inference_accepts_sparse_x_with_pre_encoder_covariates(
    pre_encoder_covariate: str,
):
    n_input = 5
    batch_size = 4
    use_continuous_covariate = pre_encoder_covariate == "continuous"
    vae = VAE(
        n_input=n_input,
        n_batch=2,
        n_hidden=8,
        n_latent=3,
        n_continuous_cov=int(use_continuous_covariate),
        encode_covariates=True,
        batch_representation=(
            "embedding" if pre_encoder_covariate == "embedding_batch" else "one-hot"
        ),
        use_batch_norm="none",
        use_layer_norm="none",
    )
    x = torch.tensor(
        [
            [1.0, 0.0, 2.0, 0.0, 1.0],
            [0.0, 3.0, 0.0, 1.0, 0.0],
            [2.0, 0.0, 1.0, 0.0, 0.0],
            [0.0, 1.0, 0.0, 2.0, 1.0],
        ]
    )
    batch_index = torch.tensor([[0], [1], [0], [1]])
    cont_covs = torch.randn(batch_size, 1) if use_continuous_covariate else None

    inference_outputs = vae._regular_inference(
        x.to_sparse_csr(),
        batch_index,
        cont_covs=cont_covs,
    )

    assert inference_outputs[MODULE_KEYS.Z_KEY].shape == (batch_size, 3)
