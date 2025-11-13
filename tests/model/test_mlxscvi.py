import sys

import pytest

from scvi.data import synthetic_iid
from scvi.model import mlxSCVI
from scvi.utils import attrdict

# the whole file should only run on macOS
pytestmark = pytest.mark.skipif(
    sys.platform != "darwin", reason="This test file runs only on macOS"
)


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

    z1 = model.get_latent_representation(give_mean=True, n_samples=1)
    assert z1.ndim == 2
    z2 = model.get_latent_representation(give_mean=False, n_samples=5)
    assert z2.ndim == 3
    assert z2.shape[0] == 5


@pytest.mark.parametrize("n_latent", [5])
@pytest.mark.parametrize("dropout_rate", [0.1])
def test_mlx_scvi_training(n_latent: int, dropout_rate: float):
    adata = synthetic_iid()
    mlxSCVI.setup_anndata(
        adata,
        batch_key="batch",
    )

    model = mlxSCVI(adata, n_latent=n_latent, dropout_rate=dropout_rate)
    assert model.module.training


@pytest.mark.parametrize("n_latent", [5])
def test_mlx_scvi_save_load(n_latent: int, save_path: str):
    adata = synthetic_iid()
    mlxSCVI.setup_anndata(
        adata,
        batch_key="batch",
    )
    model = mlxSCVI(adata, n_latent=n_latent)
    model.train(2, train_size=0.5, check_val_every_n_epoch=1)
    # z1 = model.get_latent_representation(adata)
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

    # z2 = model.get_latent_representation()
    # np.testing.assert_array_equal(z1, z2)
