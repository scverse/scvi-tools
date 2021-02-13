from scvi.data import setup_anndata, synthetic_iid
from scvi.model import SCVI
from scvi.external import SOLO


def test_solo(save_path):
    n_latent = 5
    adata = synthetic_iid(run_setup_anndata=False)
    setup_anndata(adata)
    model = SCVI(adata, n_latent=n_latent)
    model.train(1, check_val_every_n_epoch=1, train_size=0.5)

    solo = SOLO.from_scvi_model(model)
    solo.train(1)
    solo.predict()
