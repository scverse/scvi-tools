import os
from time import time

import pytest
import torch

from scvi.data import synthetic_iid
from scvi.model import CondSCVI


@pytest.mark.parametrize("n_batches", [1, 2, 3])
@pytest.mark.parametrize("encode_covariates", [True, False])
def test_condscvi_batch_key(
    save_path: str, n_batches: int, encode_covariates: bool, n_labels: int = 5
):
    adata = synthetic_iid(n_batches=n_batches, n_labels=n_labels)
    CondSCVI.setup_anndata(adata, batch_key="batch", labels_key="labels")
    model = CondSCVI(adata, encode_covariates=encode_covariates)

    model.train(max_epochs=1)
    assert model.summary_stats.n_batch == n_batches
    _ = model.get_elbo()
    _ = model.get_reconstruction_error()
    _ = model.get_latent_representation()
    _ = model.get_vamp_prior(adata)

    model.get_normalized_expression(adata)
    model.differential_expression(groupby="labels", group1="label_1")
    model.differential_expression(groupby="labels", group1="label_1", group2="label_2")

    model_path = os.path.join(save_path, __name__)
    model.save(model_path, overwrite=True, save_anndata=False)
    model = CondSCVI.load(model_path, adata=adata)


def test_condscvi_batch_key_compat_load(save_path: str):
    adata = synthetic_iid()
    model = CondSCVI.load("tests/test_data/condscvi_pre_batch", adata=adata)

    assert not hasattr(model.summary_stats, "n_batch")
    _ = model.get_latent_representation()
    _ = model.get_vamp_prior(adata)

    model_path = os.path.join(save_path, __name__)
    model.save(model_path, overwrite=True, save_anndata=False)
    model = CondSCVI.load(model_path, adata=adata)


@pytest.mark.parametrize("weight_obs", [True, False])
def test_condscvi_no_batch_key(save_path: str, weight_obs: bool):
    adata = synthetic_iid()
    CondSCVI.setup_anndata(adata, labels_key="labels")

    with pytest.raises(ValueError):
        _ = CondSCVI(adata, encode_covariates=True)

    model = CondSCVI(adata, weight_obs=weight_obs)
    model.train(max_epochs=1)
    assert not hasattr(model.summary_stats, "n_batch")
    _ = model.get_elbo()
    _ = model.get_reconstruction_error()
    _ = model.get_latent_representation()
    _ = model.get_vamp_prior(adata)

    model_path = os.path.join(save_path, __name__)
    model.save(model_path, overwrite=True, save_anndata=False)
    model = CondSCVI.load(model_path, adata=adata)


def test_cpu_gpu_condscvi():
    if torch.cuda.is_available():
        adata = synthetic_iid(10000, 500)

        CondSCVI.setup_anndata(adata, labels_key="labels")

        m = CondSCVI(adata)
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
        m2 = CondSCVI(adata)
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
