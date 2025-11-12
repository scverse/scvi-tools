from unittest import mock

import numpy as np
import pytest

from scvi.data import synthetic_iid
from scvi.model import mlxSCVI
from scvi.train import MlxTrainingPlan
from scvi.utils import attrdict


@pytest.mark.parametrize("n_latent", [5])
def test_mlx_scvi(n_latent: int):
    adata = synthetic_iid()
    mlxSCVI.setup_anndata(
        adata,
        batch_key="batch",
    )
    model = mlxSCVI(adata, n_latent=n_latent)
    model.train(2, train_size=0.5, check_val_every_n_epoch=1)
    model.get_latent_representation()

    model = mlxSCVI(adata, n_latent=n_latent, gene_likelihood="poisson")
    model.train(1, train_size=0.5)
    z1 = model.get_latent_representation(give_mean=True, n_samples=1)
    assert z1.ndim == 2
    z2 = model.get_latent_representation(give_mean=False, n_samples=15)
    assert z2.ndim == 3
    assert z2.shape[0] == 15


@pytest.mark.parametrize("n_latent", [5])
@pytest.mark.parametrize("dropout_rate", [0.1])
def test_mlx_scvi_training(n_latent: int, dropout_rate: float):
    from flax import linen as nn

    adata = synthetic_iid()
    mlxSCVI.setup_anndata(
        adata,
        batch_key="batch",
    )

    model = mlxSCVI(adata, n_latent=n_latent, dropout_rate=dropout_rate)
    assert model.module.training

    with mock.patch("scvi.module._mlxvae.nn.Dropout", wraps=nn.Dropout) as mock_dropout_cls:
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


@pytest.mark.parametrize("n_latent", [5])
def test_mlx_scvi_save_load(n_latent: int, save_path: str):
    adata = synthetic_iid()
    mlxSCVI.setup_anndata(
        adata,
        batch_key="batch",
    )
    model = mlxSCVI(adata, n_latent=n_latent)
    model.train(2, train_size=0.5, check_val_every_n_epoch=1)
    z1 = model.get_latent_representation(adata)
    model.save(save_path, overwrite=True, save_anndata=True)
    model.view_setup_args(save_path)
    model = mlxSCVI.load(save_path)
    model.get_latent_representation()

    # Load with mismatched genes.
    tmp_adata = synthetic_iid(
        n_genes=200,
    )
    with pytest.raises(ValueError):
        mlxSCVI.load(save_path, adata=tmp_adata)

    # Load with different batches.
    tmp_adata = synthetic_iid()
    tmp_adata.obs["batch"] = tmp_adata.obs["batch"].cat.rename_categories(["batch_2", "batch_3"])
    with pytest.raises(ValueError):
        mlxSCVI.load(save_path, adata=tmp_adata)

    model = mlxSCVI.load(save_path, adata=adata)
    assert "batch" in model.adata_manager.data_registry
    assert model.adata_manager.data_registry.batch == attrdict(
        {"attr_name": "obs", "attr_key": "_scvi_batch"}
    )
    assert model.is_trained is True

    z2 = model.get_latent_representation()
    np.testing.assert_array_equal(z1, z2)


def test_loss_args_mlx():
    """Test that self._loss_args is set correctly."""
    adata = synthetic_iid()
    mlxSCVI.setup_anndata(adata)
    mlx_vae = mlxSCVI(adata)
    mlx_tp = MlxTrainingPlan(mlx_vae.module)

    loss_args = [
        "tensors",
        "inference_outputs",
        "generative_outputs",
        "kl_weight",
    ]
    assert len(mlx_tp._loss_args) == len(loss_args)
    for arg in loss_args:
        assert arg in mlx_tp._loss_args
