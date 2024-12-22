from time import time

import numpy as np
import pytest
import torch

from scvi.data import synthetic_iid
from scvi.model import MULTIVI


def test_multivi():
    data = synthetic_iid()
    MULTIVI.setup_anndata(
        data,
        batch_key="batch",
    )
    vae = MULTIVI(
        data,
        n_genes=50,
        n_regions=50,
    )
    vae.train(1, save_best=False)
    vae.train(1, adversarial_mixing=False)
    vae.train(3)
    vae.get_elbo(indices=vae.validation_indices)
    vae.get_accessibility_estimates()
    vae.get_accessibility_estimates(normalize_cells=True)
    vae.get_accessibility_estimates(normalize_regions=True)
    vae.get_normalized_expression()
    vae.get_library_size_factors()
    vae.get_region_factors()
    vae.get_reconstruction_error(indices=vae.validation_indices)
    vae.get_latent_representation()
    vae.differential_accessibility(groupby="labels", group1="label_1")
    vae.differential_expression(groupby="labels", group1="label_1")

    # Test with size factor
    data = synthetic_iid()
    data.obs["size_factor"] = np.random.randint(1, 5, size=(data.shape[0],))
    MULTIVI.setup_anndata(data, batch_key="batch", size_factor_key="size_factor")
    vae = MULTIVI(
        data,
        n_genes=50,
        n_regions=50,
    )
    vae.train(3)

    # Test with modality weights and penalties
    data = synthetic_iid()
    MULTIVI.setup_anndata(data, batch_key="batch")
    vae = MULTIVI(data, n_genes=50, n_regions=50, modality_weights="cell")
    vae.train(3)
    vae = MULTIVI(data, n_genes=50, n_regions=50, modality_weights="universal")
    vae.train(3)
    vae = MULTIVI(data, n_genes=50, n_regions=50, modality_penalty="MMD")
    vae.train(3)

    # Test with non-zero protein data
    data = synthetic_iid()
    MULTIVI.setup_anndata(
        data,
        batch_key="batch",
        protein_expression_obsm_key="protein_expression",
        protein_names_uns_key="protein_names",
    )
    vae = MULTIVI(
        data,
        n_genes=50,
        n_regions=50,
        modality_weights="cell",
    )
    assert vae.n_proteins == data.obsm["protein_expression"].shape[1]
    vae.train(3)


def test_multivi_single_batch():
    data = synthetic_iid(n_batches=1)
    MULTIVI.setup_anndata(data, batch_key="batch")
    vae = MULTIVI(data, n_genes=50, n_regions=50)
    with pytest.warns(UserWarning):
        vae.train(3)


@pytest.mark.optional
def test_cpu_gpu_multivi():
    if torch.cuda.is_available():
        adata = synthetic_iid(10000, 500)

        MULTIVI.setup_anndata(
            adata,
            batch_key="batch",
        )

        m = MULTIVI(adata, n_genes=50, n_regions=50)
        training_start_time = time()
        m.train(
            accelerator="cpu",
            batch_size=5000,
            max_epochs=100,
            train_size=0.9,
            plan_kwargs={"n_epochs_kl_warmup": 100, "compile": False},
            datasplitter_kwargs={"drop_last": True},
        )
        print(f"CPU Training finished, took {time() - training_start_time:.2f}s")
        m.get_latent_representation()
        m.get_elbo()
        m.get_reconstruction_error()

        # run the exact same thing on GPU:
        m2 = MULTIVI(adata, n_genes=50, n_regions=50)
        training_start_time2 = time()
        m2.train(
            accelerator="cuda",
            batch_size=5000,
            max_epochs=100,
            train_size=0.9,
            plan_kwargs={"n_epochs_kl_warmup": 100, "compile": True},
            datasplitter_kwargs={"drop_last": True},
        )
        print(f"Compile + GPU Training finished, took {time() - training_start_time2:.2f}s")
        m2.get_latent_representation()
        m2.get_elbo()
        m2.get_reconstruction_error()
