import torch

import scvi
from scvi.data import synthetic_iid
from scvi.model import SCVI
from scvi.train._callbacks import SaveBestState


def test_save_best_state_callback(save_path):
    n_latent = 5
    adata = synthetic_iid()
    SCVI.setup_anndata(adata, batch_key="batch", labels_key="labels")
    model = SCVI(adata, n_latent=n_latent)
    callbacks = [SaveBestState(verbose=True)]
    model.train(3, check_val_every_n_epoch=1, train_size=0.5, callbacks=callbacks)


def test_set_seed(save_path):
    scvi.settings.seed = 1
    n_latent = 5
    adata = synthetic_iid()
    SCVI.setup_anndata(adata, batch_key="batch", labels_key="labels")
    model1 = SCVI(adata, n_latent=n_latent)
    model1.train(1)
    scvi.settings.seed = 1
    model2 = SCVI(adata, n_latent=n_latent)
    model2.train(1)
    assert torch.equal(
        model1.module.z_encoder.encoder.fc_layers[0][0].weight,
        model2.module.z_encoder.encoder.fc_layers[0][0].weight,
    )
