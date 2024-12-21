import numpy as np
import pytest

from scvi.data import synthetic_iid
from scvi.external import Decipher
from scvi.external.decipher.utils import Trajectory, rotate_decipher_components


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


def test_decipher_trajectory(adata):
    Decipher.setup_anndata(adata)
    model = Decipher(adata)
    model.train(
        max_epochs=2,
        check_val_every_n_epoch=1,
        train_size=0.5,
        early_stopping=True,
    )

    v = model.get_latent_representation(give_z=False)
    adata.obsm["decipher_v"] = v

    trajectory = Trajectory.from_dict(
        {
            "rep_key": "decipher_v",
            "cluster_locations": v[
                np.random.choice(np.arange(v.shape[0]), size=10, replace=False), :
            ],
            "cluster_ids": np.arange(10),
            "density": 50,
        }
    )

    adata.obs["cluster_ids"] = np.random.randint(0, 10, size=adata.n_obs)
    decipher_time = model.compute_decipher_time(
        adata,
        cluster_obs_key="cluster_ids",
        trajectory=trajectory,
        n_neighbors=10,
    )
    assert decipher_time.shape == (adata.n_obs,)

    gene_patterns = model.compute_gene_patterns(
        adata,
        trajectory=trajectory,
        l_scale=10_000,
        n_samples=100,
    )
    assert gene_patterns["mean"].shape[1] == adata.n_vars
    assert gene_patterns["q25"].shape[1] == adata.n_vars
    assert gene_patterns["q75"].shape[1] == adata.n_vars
    n_times = gene_patterns["mean"].shape[0]
    assert gene_patterns["times"].shape == (n_times,)
