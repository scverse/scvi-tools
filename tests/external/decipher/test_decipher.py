import pytest
import numpy as np

from scvi.data import synthetic_iid
from scvi.external import Decipher
from scvi.external.decipher.utils import rotate_decipher_components


@pytest.fixture(scope="session")
def adata():
    adata = synthetic_iid()
    return adata


def test_decipher_train(adata):
    Decipher.setup_anndata(adata)
    model = Decipher(adata)
    model.train(
        max_epochs=2,
        check_val_every_n_epoch=1,
        train_size=0.5,
        early_stopping=True,
    )


def test_decipher_model_methods(adata):
    Decipher.setup_anndata(adata)
    model = Decipher(adata)
    model.train(
        max_epochs=2,
        check_val_every_n_epoch=1,
        train_size=0.5,
        early_stopping=True,
    )
    v = model.get_latent_representation(give_z=False)
    z = model.get_latent_representation(give_z=True)
    assert v.shape == (adata.n_obs, model.module.dim_v)
    assert z.shape == (adata.n_obs, model.module.dim_z)

    v_obsm_key = "decipher_v"
    z_obsm_key = "decipher_z"
    adata.obsm[v_obsm_key] = v
    adata.obsm[z_obsm_key] = z

    v1_obs_col = "batch"
    v2_obs_col = "label"
    rotated_v, rotated_z, rotation = rotate_decipher_components(
        adata,
        v_obsm_key=v_obsm_key,
        z_obsm_key=z_obsm_key,
        v1_obs_col=v1_obs_col,
        v1_order=["batch_0", "batch_1"],
        auto_flip_decipher_z=True,
    )
    assert rotated_v.shape == (adata.n_obs, model.module.dim_v)
    assert not np.allclose(rotated_v, v)
    assert rotated_z.shape == (adata.n_obs, model.module.dim_z)
    assert not np.allclose(rotated_z, z)
    assert rotation.shape == (model.module.dim_v, model.module.dim_v)

    imputed_expr, v_cov, z_cov = model.compute_imputed_gene_expression(
        adata,
        compute_covariances=True,
        v_obsm_key=v_obsm_key,
        z_obsm_key=z_obsm_key,
    )
    assert imputed_expr.shape == (adata.n_obs, adata.n_vars)
    assert v_cov.shape == (adata.n_vars, model.module.dim_v)
    assert z_cov.shape == (adata.n_vars, model.module.dim_z)
