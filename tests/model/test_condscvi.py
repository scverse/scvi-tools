import pytest

from scvi.data import synthetic_iid
from scvi.model import CondSCVI


def test_condscvi():
    dataset = synthetic_iid(
        n_labels=5,
    )
    CondSCVI.setup_anndata(
        dataset,
        labels_key="labels",
    )
    model = CondSCVI(dataset)
    model.train(1, train_size=1)
    model.get_latent_representation()
    model.get_vamp_prior(dataset)

    model = CondSCVI(dataset, weight_obs=True)
    model.train(1, train_size=1)
    model.get_latent_representation()
    model.get_vamp_prior(dataset)


def test_condscvi_no_batch_key(n_batches: int = 3, n_labels: int = 5):
    adata = synthetic_iid(n_batches=n_batches, n_labels=n_labels)
    CondSCVI.setup_anndata(adata, labels_key="labels")

    with pytest.raises(ValueError):
        _ = CondSCVI(adata, encode_covariates=True)

    model = CondSCVI(adata)

    model.train(max_epochs=1)


def test_condscvi_batch_key(n_batches: int = 3, n_labels: int = 5):
    adata = synthetic_iid(n_batches=n_batches, n_labels=n_labels)
    CondSCVI.setup_anndata(adata, batch_key="batch", labels_key="labels")
    model = CondSCVI(adata)

    assert model.summary_stats.n_batch == n_batches
