import sys
import tempfile

import numpy as np
import pytest

from scvi.data import synthetic_iid
from scvi.utils import attrdict

# the whole file should only run on macOS
pytestmark = pytest.mark.skipif(
    sys.platform != "darwin", reason="This test file runs only on macOS"
)


@pytest.mark.parametrize("n_latent", [5])
def test_mlx_scvi(n_latent: int):
    from scvi.model import mlxSCVI

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
    from scvi.model import mlxSCVI

    adata = synthetic_iid()
    mlxSCVI.setup_anndata(
        adata,
        batch_key="batch",
    )

    model = mlxSCVI(adata, n_latent=n_latent, dropout_rate=dropout_rate)
    assert model.module.training


@pytest.mark.parametrize("n_latent", [5])
def test_mlx_scvi_save_load(n_latent: int, save_path: str):
    from scvi.model import mlxSCVI

    adata = synthetic_iid()
    mlxSCVI.setup_anndata(
        adata,
        batch_key="batch",
    )
    model = mlxSCVI(adata, n_latent=n_latent)
    model.train(2, train_size=0.5, check_val_every_n_epoch=1)
    z1 = model.get_latent_representation(adata)

    # Verify trained weights are MLX arrays with expected shapes.
    # MlxDense stores weight as (in_features, out_features).
    import mlx.core as mx

    enc_w = model.module.encoder.dense1.weight
    assert isinstance(enc_w, mx.array), "Encoder weights should be MLX arrays after training"
    assert enc_w.shape == (adata.n_vars, 128), (
        f"Encoder dense1 weight: expected ({adata.n_vars}, 128), got {enc_w.shape}"
    )
    dec_w = model.module.decoder.dense1.weight
    assert isinstance(dec_w, mx.array), "Decoder weights should be MLX arrays after training"
    assert dec_w.shape == (n_latent, 128), (
        f"Decoder dense1 weight: expected ({n_latent}, 128), got {dec_w.shape}"
    )

    model.save(save_path, overwrite=True, save_anndata=True, prefix="mlx_")
    model.view_setup_args(save_path)
    model = mlxSCVI.load(save_path, prefix="mlx_")
    model.get_latent_representation()

    # Load with mismatched genes.
    tmp_adata = synthetic_iid(
        n_genes=200,
    )
    with pytest.raises(ValueError):
        mlxSCVI.load(save_path, adata=tmp_adata, prefix="mlx_")

    # Load with different batches.
    tmp_adata = synthetic_iid()
    tmp_adata.obs["batch"] = tmp_adata.obs["batch"].cat.rename_categories(["batch_2", "batch_3"])
    with pytest.raises(ValueError):
        mlxSCVI.load(save_path, adata=tmp_adata, prefix="mlx_")

    model = mlxSCVI.load(save_path, adata=adata, prefix="mlx_")
    assert "batch" in model.adata_manager.data_registry
    assert model.adata_manager.data_registry.batch == attrdict(
        {"attr_name": "obs", "attr_key": "_scvi_batch"}
    )
    assert model.is_trained is True

    z2 = model.get_latent_representation()
    np.testing.assert_array_equal(z1, z2)

    # After loading, weights must still be MLX arrays (not plain lists or tensors).
    loaded_enc_w = model.module.encoder.dense1.weight
    assert isinstance(loaded_enc_w, mx.array), (
        "Encoder weights should remain MLX arrays after save/load"
    )
    loaded_dec_w = model.module.decoder.dense1.weight
    assert isinstance(loaded_dec_w, mx.array), (
        "Decoder weights should remain MLX arrays after save/load"
    )

    # Trying to load an mlxSCVI-saved model using the plain SCVI class must fail.
    from scvi.model import SCVI

    with pytest.raises(ValueError):
        SCVI.load(save_path, adata=adata, prefix="mlx_")

    # Trying to load a plain SCVI model using mlxSCVI must also fail.
    scvi_adata = synthetic_iid()
    SCVI.setup_anndata(scvi_adata, batch_key="batch")
    scvi_model = SCVI(scvi_adata, n_latent=n_latent)
    scvi_model.train(2, train_size=0.5)
    with tempfile.TemporaryDirectory() as scvi_tmp_path:
        scvi_model.save(scvi_tmp_path, overwrite=True)
        with pytest.raises(ValueError):
            mlxSCVI.load(scvi_tmp_path, adata=scvi_adata)

    # Now load new data/model with the old model and make training step
    adata_new = synthetic_iid()
    mlxSCVI.setup_anndata(
        adata_new,
        batch_key="batch",
    )
    model_new = mlxSCVI.load(save_path, adata=adata_new, prefix="mlx_")
    assert "batch" in model_new.adata_manager.data_registry
    assert model_new.adata_manager.data_registry.batch == attrdict(
        {"attr_name": "obs", "attr_key": "_scvi_batch"}
    )
    assert model_new.is_trained is True
    model_new.train(2, train_size=0.5, check_val_every_n_epoch=1)


@pytest.mark.parametrize("n_latent", [5])
def test_mlx_scvi_loss_decreases(n_latent: int):
    """Verify loss actually decreases during training."""
    import mlx.core as mx

    from scvi.model import mlxSCVI

    adata = synthetic_iid()
    mlxSCVI.setup_anndata(adata, batch_key="batch")
    model = mlxSCVI(adata, n_latent=n_latent)

    # Compute loss before training
    scdl = model._make_data_loader(adata=adata, batch_size=128, iter_ndarray=True)
    batch = next(iter(scdl))
    tensors = {k: mx.array(v) for k, v in batch.items()}
    model.module.train()
    _, _, loss_before = model.module(tensors)
    early_loss = float(loss_before.loss)

    # Train for several epochs
    model.train(50, train_size=1.0)

    # Compute loss after training
    model.module.eval()
    _, _, loss_after = model.module(tensors)
    late_loss = float(loss_after.loss)

    assert early_loss > 0, f"Loss should be positive, got {early_loss:.4f}"
    assert late_loss > 0, f"Loss should be positive, got {late_loss:.4f}"
    assert late_loss < early_loss, (
        f"Loss did not decrease: early={early_loss:.4f}, late={late_loss:.4f}"
    )


@pytest.mark.parametrize("n_latent", [5])
def test_mlx_scvi_poisson(n_latent: int):
    """Test that Poisson likelihood trains without errors."""
    import mlx.core as mx

    from scvi.model import mlxSCVI

    adata = synthetic_iid()
    mlxSCVI.setup_anndata(adata, batch_key="batch")
    model = mlxSCVI(adata, n_latent=n_latent, gene_likelihood="poisson")
    model.train(5, train_size=1.0)
    z = model.get_latent_representation()
    assert z.shape == (adata.n_obs, n_latent)

    # Verify loss is positive
    scdl = model._make_data_loader(adata=adata, batch_size=128, iter_ndarray=True)
    batch = next(iter(scdl))
    tensors = {k: mx.array(v) for k, v in batch.items()}
    model.module.eval()
    _, _, loss_output = model.module(tensors)
    assert float(loss_output.loss) > 0, "Poisson loss should be positive"


def test_multiple_covariates_mlxscvi():
    """Test that mlxSCVI can handle multiple categorical and continuous covariates."""
    from scvi.model import mlxSCVI

    adata = synthetic_iid()
    adata.obs["cont1"] = np.random.normal(size=(adata.shape[0],))
    adata.obs["cont2"] = np.random.normal(size=(adata.shape[0],))
    adata.obs["cat1"] = np.random.randint(0, 5, size=(adata.shape[0],))
    adata.obs["cat2"] = np.random.randint(0, 5, size=(adata.shape[0],))

    mlxSCVI.setup_anndata(
        adata,
        batch_key="batch",
        labels_key="labels",
        continuous_covariate_keys=["cont1", "cont2"],
        categorical_covariate_keys=["cat1", "cat2"],
    )
    m = mlxSCVI(adata)
    m.train(1)
    m.get_latent_representation()
