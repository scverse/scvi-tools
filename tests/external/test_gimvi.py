import numpy as np
import pytest

import scvi
from scvi.data import synthetic_iid
from scvi.external import GIMVI


def test_saving_and_loading(save_path):
    prefix = "GIMVI_"
    adata = synthetic_iid()
    GIMVI.setup_anndata(
        adata,
        batch_key="batch",
    )
    adata2 = synthetic_iid()
    GIMVI.setup_anndata(
        adata2,
        batch_key="batch",
    )

    # GIMVI
    model = GIMVI(adata, adata2)
    model.train(3, train_size=0.5)
    z1 = model.get_latent_representation([adata])
    z2 = model.get_latent_representation([adata])
    np.testing.assert_array_equal(z1, z2)
    model.save(save_path, overwrite=True, save_anndata=True, prefix=prefix)
    model = GIMVI.load(save_path, prefix=prefix)
    model.get_latent_representation()
    tmp_adata = scvi.data.synthetic_iid(n_genes=200)
    tmp_adata2 = scvi.data.synthetic_iid(n_genes=200)
    with pytest.raises(ValueError):
        GIMVI.load(
            save_path, adata_seq=tmp_adata, adata_spatial=tmp_adata2, prefix=prefix
        )
    model = GIMVI.load(save_path, adata_seq=adata, adata_spatial=adata2, prefix=prefix)
    z2 = model.get_latent_representation([adata])
    np.testing.assert_array_equal(z1, z2)
    model = GIMVI.load(
        save_path,
        adata_seq=adata,
        adata_spatial=adata2,
        use_gpu=False,
        prefix=prefix,
    )
    z2 = model.get_latent_representation([adata])
    np.testing.assert_almost_equal(z1, z2, decimal=3)
    assert model.is_trained is True


def test_gimvi():
    adata_seq = synthetic_iid()
    adata_spatial = synthetic_iid()
    GIMVI.setup_anndata(
        adata_seq,
        batch_key="batch",
        labels_key="labels",
    )
    GIMVI.setup_anndata(
        adata_spatial,
        batch_key="batch",
        labels_key="labels",
    )
    model = GIMVI(adata_seq, adata_spatial, n_latent=10)
    assert hasattr(model.module, "library_log_means_0") and not hasattr(
        model.module, "library_log_means_1"
    )
    model.train(1, check_val_every_n_epoch=1, train_size=0.5)
    model.get_latent_representation()
    model.get_imputed_values()

    adata_spatial.var_names += "asdf"
    GIMVI.setup_anndata(
        adata_spatial,
        batch_key="batch",
        labels_key="labels",
    )
    with pytest.raises(ValueError):
        model = GIMVI(adata_seq, adata_spatial)


def test_gimvi_model_library_size():
    adata_seq = synthetic_iid()
    adata_spatial = synthetic_iid()
    GIMVI.setup_anndata(
        adata_seq,
        batch_key="batch",
        labels_key="labels",
    )
    GIMVI.setup_anndata(
        adata_spatial,
        batch_key="batch",
        labels_key="labels",
    )
    model = GIMVI(
        adata_seq, adata_spatial, model_library_size=[True, True], n_latent=10
    )
    assert hasattr(model.module, "library_log_means_0") and hasattr(
        model.module, "library_log_means_1"
    )
    model.train(1, check_val_every_n_epoch=1, train_size=0.5)
    model.get_latent_representation()
    model.get_imputed_values()
