from unittest import mock

import numpy as np
import pytest
from flax import linen as nn

from scvi.data import synthetic_iid
from scvi.model import JaxSCVI
from scvi.utils import attrdict


def test_jax_scvi(n_latent=5):
    adata = synthetic_iid()
    JaxSCVI.setup_anndata(
        adata,
        batch_key="batch",
    )
    model = JaxSCVI(adata, n_latent=n_latent)
    model.train(2, train_size=0.5, check_val_every_n_epoch=1)
    model.get_latent_representation()

    model = JaxSCVI(adata, n_latent=n_latent, gene_likelihood="poisson")
    model.train(1, train_size=0.5)
    z1 = model.get_latent_representation(give_mean=True, n_samples=1)
    assert z1.ndim == 2
    z2 = model.get_latent_representation(give_mean=False, n_samples=15)
    assert (z2.ndim == 3) and (z2.shape[0] == 15)


def test_jax_scvi_training(n_latent: int = 5, dropout_rate: float = 0.1):
    adata = synthetic_iid()
    JaxSCVI.setup_anndata(
        adata,
        batch_key="batch",
    )

    model = JaxSCVI(adata, n_latent=n_latent, dropout_rate=dropout_rate)
    assert model.module.training

    with mock.patch(
        "scvi.module._jaxvae.nn.Dropout", wraps=nn.Dropout
    ) as mock_dropout_cls:
        mock_dropout = mock.Mock()
        mock_dropout.side_effect = lambda h, **_kwargs: h
        mock_dropout_cls.return_value = mock_dropout
        model.train(1, train_size=0.5, check_val_every_n_epoch=1)

        assert not model.module.training
        mock_dropout_cls.assert_called()
        mock_dropout.assert_has_calls(
            12 * [mock.call(mock.ANY, deterministic=False)]
            + 8 * [mock.call(mock.ANY, deterministic=True)]
        )


def test_jax_scvi_save_load(save_path: str, n_latent: int = 5):
    adata = synthetic_iid()
    JaxSCVI.setup_anndata(
        adata,
        batch_key="batch",
    )
    model = JaxSCVI(adata, n_latent=n_latent)
    model.train(2, train_size=0.5, check_val_every_n_epoch=1)
    z1 = model.get_latent_representation(adata)
    model.save(save_path, overwrite=True, save_anndata=True)
    model.view_setup_args(save_path)
    model = JaxSCVI.load(save_path)
    model.get_latent_representation()

    # Load with mismatched genes.
    tmp_adata = synthetic_iid(
        n_genes=200,
    )
    with pytest.raises(ValueError):
        JaxSCVI.load(save_path, adata=tmp_adata)

    # Load with different batches.
    tmp_adata = synthetic_iid()
    tmp_adata.obs["batch"] = tmp_adata.obs["batch"].cat.rename_categories(
        ["batch_2", "batch_3"]
    )
    with pytest.raises(ValueError):
        JaxSCVI.load(save_path, adata=tmp_adata)

    model = JaxSCVI.load(save_path, adata=adata)
    assert "batch" in model.adata_manager.data_registry
    assert model.adata_manager.data_registry.batch == attrdict(
        {"attr_name": "obs", "attr_key": "_scvi_batch"}
    )
    assert model.is_trained is True

    z2 = model.get_latent_representation()
    np.testing.assert_array_equal(z1, z2)
