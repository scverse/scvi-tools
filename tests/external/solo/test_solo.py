import numpy as np
import pytest

from scvi.data import synthetic_iid
from scvi.external import SOLO
from scvi.model import SCVI


@pytest.mark.parametrize("soft", [True, False])
@pytest.mark.parametrize("return_logits", [True, False])
def test_solo(soft: bool, return_logits: bool):
    n_latent = 5
    adata = synthetic_iid()
    SCVI.setup_anndata(adata)
    model = SCVI(adata, n_latent=n_latent)
    model.train(1, check_val_every_n_epoch=1, train_size=0.5)

    solo = SOLO.from_scvi_model(model)
    solo.train(1, check_val_every_n_epoch=1, train_size=0.9)
    assert "validation_loss" in solo.history.keys()
    predictions = solo.predict(soft=soft, return_logits=return_logits)
    if soft:
        preds = predictions.to_numpy()
        assert preds.shape == (adata.n_obs, 2)
        if not return_logits:
            assert np.allclose(preds.sum(axis=-1), 1)
        else:
            assert not np.allclose(preds.sum(axis=-1), 1)

    bdata = synthetic_iid()
    solo = SOLO.from_scvi_model(model, bdata)
    solo.train(1, check_val_every_n_epoch=1, train_size=0.9)
    assert "validation_loss" in solo.history.keys()
    predictions = solo.predict(soft=soft)
    if soft:
        preds = predictions.to_numpy()
        assert preds.shape == (adata.n_obs, 2)
        assert np.allclose(preds.sum(axis=-1), 1)


def test_solo_multiple_batch():
    n_latent = 5
    adata = synthetic_iid()
    adata.layers["my_layer"] = adata.X.copy()
    SCVI.setup_anndata(adata, layer="my_layer", batch_key="batch")
    model = SCVI(adata, n_latent=n_latent)
    model.train(1, check_val_every_n_epoch=1, train_size=0.5)

    solo = SOLO.from_scvi_model(model, restrict_to_batch="batch_0")
    solo.train(1, check_val_every_n_epoch=1, train_size=0.9)
    assert "validation_loss" in solo.history.keys()
    solo.predict()


def test_solo_scvi_labels():
    n_latent = 5
    adata = synthetic_iid()
    SCVI.setup_anndata(adata, labels_key="labels")
    model = SCVI(adata, n_latent=n_latent)
    model.train(1, check_val_every_n_epoch=1, train_size=0.5)

    solo = SOLO.from_scvi_model(model)
    solo.train(1, check_val_every_n_epoch=1, train_size=0.9)
    assert "validation_loss" in solo.history.keys()
    solo.predict()


def test_solo_from_scvi_errors():
    adata = synthetic_iid()
    adata.obs["continuous_covariate"] = np.random.normal(size=(adata.n_obs, 1))
    adata.obs["categorical_covariate"] = np.random.choice(["a", "b", "c"], size=(adata.n_obs, 1))

    # no batch key, restrict_to_batch
    SCVI.setup_anndata(adata, labels_key="labels")
    model = SCVI(adata)
    model.train(max_epochs=1)
    with pytest.raises(ValueError):
        _ = SOLO.from_scvi_model(model, restrict_to_batch="batch_0")

    # continuous covariate
    SCVI.setup_anndata(
        adata,
        batch_key="batch",
        labels_key="labels",
        continuous_covariate_keys=["continuous_covariate"],
    )
    model = SCVI(adata)
    model.train(max_epochs=1)

    with pytest.raises(ValueError):
        _ = SOLO.from_scvi_model(model, restrict_to_batch="batch_0")

    # categorical covariate
    SCVI.setup_anndata(
        adata,
        batch_key="batch",
        labels_key="labels",
        categorical_covariate_keys=["categorical_covariate"],
    )
    model = SCVI(adata)
    model.train(max_epochs=1)

    with pytest.raises(ValueError):
        _ = SOLO.from_scvi_model(model, restrict_to_batch="batch_0")
