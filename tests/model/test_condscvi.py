import os

import numpy as np
import pytest

from scvi.data import synthetic_iid
from scvi.model import CondSCVI


@pytest.mark.parametrize("n_batches", [1, 2, 3])
@pytest.mark.parametrize("encode_covariates", [True, False])
@pytest.mark.parametrize("n_labels", [5])
def test_condscvi_batch_key(
    n_batches: int, encode_covariates: bool, n_labels: int, save_path: str
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


def test_condscvi_fine_celltype(save_path: str):
    adata = synthetic_iid(n_batches=5, n_labels=5)
    adata.obs["fine_labels"] = [i + str(np.random.randint(2)) for i in adata.obs["labels"]]
    CondSCVI.setup_anndata(
        adata, batch_key="batch", labels_key="labels", fine_labels_key="fine_labels"
    )
    model = CondSCVI(adata, encode_covariates=True)

    model.train(max_epochs=1)
    model.predict()
    model.predict(adata=adata)
    model.predict(adata, soft=True, use_posterior_mean=False)


def test_condscvi_batch_key_compat_load(save_path: str):
    adata = synthetic_iid(n_batches=1, n_labels=5)
    model = CondSCVI.load("tests/test_data/condscvi_pre_batch", adata=adata)

    # assert not hasattr(model.summary_stats, "n_batch")
    _ = model.get_latent_representation()
    _ = model.get_vamp_prior(adata)

    model_path = os.path.join(save_path, __name__)
    model.save(model_path, overwrite=True, save_anndata=False)
    model = CondSCVI.load(model_path, adata=adata)


@pytest.mark.parametrize("weight_obs", [True, False])
def test_condscvi_no_batch_key(weight_obs: bool, save_path: str):
    adata = synthetic_iid()
    CondSCVI.setup_anndata(adata, labels_key="labels")

    # with pytest.raises(ValueError):
    #     _ = CondSCVI(adata, encode_covariates=True)

    model = CondSCVI(adata, weight_obs=weight_obs)
    model.train(max_epochs=1)
    # assert not hasattr(model.summary_stats, "n_batch")
    _ = model.get_elbo()
    _ = model.get_reconstruction_error()
    _ = model.get_latent_representation()
    _ = model.get_vamp_prior(adata)

    model_path = os.path.join(save_path, __name__)
    model.save(model_path, overwrite=True, save_anndata=False)
    model = CondSCVI.load(model_path, adata=adata)
