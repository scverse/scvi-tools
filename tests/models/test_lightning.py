from scvi.data import synthetic_iid
from scvi.model import SCVI
from scvi.lightning._callbacks import SaveBestState


def test_save_best_state_callback(save_path):

    n_latent = 5
    adata = synthetic_iid()
    model = SCVI(adata, n_latent=n_latent)
    callbacks = [SaveBestState(verbose=True)]
    model.train(3, check_val_every_n_epoch=1, train_size=0.5, callbacks=callbacks)
