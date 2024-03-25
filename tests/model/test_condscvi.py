import pytest

from scvi.data import synthetic_iid
from scvi.model import CondSCVI


@pytest.mark.parametrize("n_batches", [1, 2, 3])
def test_condscvi_batch_key(n_batches: int, n_labels: int = 5):
    adata = synthetic_iid(n_batches=n_batches, n_labels=n_labels)
    CondSCVI.setup_anndata(adata, batch_key="batch", labels_key="labels")
    model = CondSCVI(adata)

    model.train(max_epochs=1)
    assert model.summary_stats.n_batch == n_batches
    _ = model.get_latent_representation()
    _ = model.get_vamp_prior(adata)


@pytest.mark.parametrize("weight_obs", [True, False])
def test_condscvi_no_batch_key(weight_obs: bool):
    adata = synthetic_iid()
    CondSCVI.setup_anndata(adata, labels_key="labels")

    with pytest.raises(ValueError):
        _ = CondSCVI(adata, encode_covariates=True)

    model = CondSCVI(adata, weight_obs=weight_obs)
    model.train(max_epochs=1)
    assert not hasattr(model.summary_stats, "n_batch")
    _ = model.get_latent_representation()
    _ = model.get_vamp_prior(adata)
