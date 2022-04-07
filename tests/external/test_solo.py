from scvi.data import synthetic_iid
from scvi.external import SOLO
from scvi.model import SCVI


def test_solo():
    n_latent = 5
    adata = synthetic_iid()
    SCVI.setup_anndata(adata)
    model = SCVI(adata, n_latent=n_latent)
    model.train(1, check_val_every_n_epoch=1, train_size=0.5)

    solo = SOLO.from_scvi_model(model)
    solo.train(1, check_val_every_n_epoch=1, train_size=0.9)
    assert "validation_loss" in solo.history.keys()
    solo.predict()

    bdata = synthetic_iid()
    solo = SOLO.from_scvi_model(model, bdata)
    solo.train(1, check_val_every_n_epoch=1, train_size=0.9)
    assert "validation_loss" in solo.history.keys()
    solo.predict()


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
