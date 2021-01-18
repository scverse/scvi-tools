import pandas as pd
import numpy as np
import pytest

import scvi
from scvi.data import synthetic_iid
from scvi.external import GIMVI


def test_saving_and_loading(save_path):
    adata = synthetic_iid()

    # GIMVI
    model = GIMVI(adata, adata)
    model.train(1, train_size=0.5)
    z1 = model.get_latent_representation([adata])
    z2 = model.get_latent_representation([adata])
    np.testing.assert_array_equal(z1, z2)
    model.save(save_path, overwrite=True, save_anndata=True)
    model = GIMVI.load(save_path)
    model.get_latent_representation()
    tmp_adata = scvi.data.synthetic_iid(n_genes=200)
    with pytest.raises(ValueError):
        GIMVI.load(save_path, tmp_adata, tmp_adata)
    model = GIMVI.load(save_path, adata, adata)
    z2 = model.get_latent_representation([adata])
    np.testing.assert_array_equal(z1, z2)
    assert model.is_trained is True


def test_gimvi():
    adata_seq = synthetic_iid()
    adata_spatial = synthetic_iid()
    model = GIMVI(adata_seq, adata_spatial, n_latent=10)
    model.train(1, check_val_every_n_epoch=1, train_size=0.5)
    model.get_latent_representation()
    model.get_imputed_values()

    assert len(model.history["elbo_train_0"]) == 2
    assert len(model.history["elbo_train_1"]) == 2
    assert len(model.history["elbo_test_0"]) == 2
    assert len(model.history["elbo_test_1"]) == 2

    trainer = model.trainer
    results = pd.DataFrame(
        trainer.get_loss_magnitude(),
        index=["reconstruction", "kl_divergence", "discriminator"],
        columns=["Sequencing", "Spatial"],
    )
    results.columns.name = "Dataset"
    results.index.name = "Loss"
    trainer.get_discriminator_confusion()
    adata_spatial.var_names += "asdf"
    with pytest.raises(ValueError):
        model = GIMVI(adata_seq, adata_spatial)
