from scvi.data import setup_anndata, synthetic_iid
from scvi.external import SOLO
from scvi.model import SCVI


def test_solo(save_path):
    n_latent = 5
    adata = synthetic_iid(run_setup_anndata=False)
    setup_anndata(adata)
    model = SCVI(adata, n_latent=n_latent)
    model.train(1, check_val_every_n_epoch=1, train_size=0.5)

    solo = SOLO.from_scvi_model(model)
    solo.train(1, check_val_every_n_epoch=1, train_size=0.9)
    assert "validation_loss" in solo.history.keys()
    solo.predict()
