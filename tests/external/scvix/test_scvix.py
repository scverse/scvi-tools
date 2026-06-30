import os

import numpy as np
import pytest

import scvi
from scvi.data import synthetic_iid
from scvi.data._constants import ADATA_MINIFY_TYPE
from scvi.data._utils import _is_minified
from scvi.external import SCVIX
from scvi.model.base import BaseMinifiedModeModelClass


def assert_approx_equal(a, b):
    # Allclose because on GPU, the values are not exactly the same
    # as some values are moved to cpu during data minification
    np.testing.assert_allclose(a, b, rtol=3e-1, atol=5e-1)


@pytest.mark.parametrize("prior", ["gaussian", "mog", "vamp", "mog_celltype"])
def test_scvix(prior: str):
    adata = synthetic_iid(batch_size=100)
    SCVIX.setup_anndata(adata, batch_key="batch", assay_key="batch", labels_key="labels")
    model = SCVIX(adata, prior=prior)
    model.train(max_epochs=1)
    model.get_latent_representation()
    model.get_normalized_expression()
    model.get_normalized_expression(transform_batch="batch_1")
    model.get_normalized_expression(n_samples=2)
    model.get_elbo(indices=model.validation_indices)
    model.get_reconstruction_error(indices=model.validation_indices)
    model.differential_expression(groupby="labels", group1="label_1")


@pytest.mark.parametrize("dispersion", ["gene", "gene-batch", "gene-assay", "gene-cell"])
def test_scvix_dispersion(dispersion: str):
    adata = synthetic_iid(batch_size=100)
    SCVIX.setup_anndata(adata, batch_key="batch", assay_key="batch", labels_key="labels")
    model = SCVIX(adata, dispersion=dispersion)
    model.train(max_epochs=1)
    model.get_normalized_expression()


def test_scvix_transform_assay():
    adata = synthetic_iid(batch_size=100)
    # Use batch as both batch and assay (two categories: batch_0, batch_1)
    SCVIX.setup_anndata(adata, batch_key="batch", assay_key="batch")
    model = SCVIX(adata)
    model.train(max_epochs=1)

    # baseline: no transformation
    expr_default = model.get_normalized_expression(return_numpy=True)

    # transform by integer index (assay 0)
    scvi.settings.seed = 42
    expr_int = model.get_normalized_expression(transform_assay=0, return_numpy=True)
    assert expr_int.shape == expr_default.shape

    # transform by string name — must match integer result (same seed, same code)
    scvi.settings.seed = 42
    expr_str = model.get_normalized_expression(transform_assay="batch_0", return_numpy=True)
    np.testing.assert_array_equal(expr_int, expr_str)

    # transform_assay and transform_batch can be combined
    expr_both = model.get_normalized_expression(
        transform_assay="batch_0", transform_batch="batch_1", return_numpy=True
    )
    assert expr_both.shape == expr_default.shape

    # different assay indices produce different outputs
    expr_assay1 = model.get_normalized_expression(transform_assay=1, return_numpy=True)
    assert not np.allclose(expr_int, expr_assay1)

    # invalid assay name raises ValueError
    with pytest.raises(ValueError, match="not a valid assay category"):
        model.get_normalized_expression(transform_assay="nonexistent_assay")


def test_scvix_encode_covariates():
    adata = synthetic_iid(batch_size=100)
    SCVIX.setup_anndata(adata, batch_key="batch", assay_key="batch", labels_key="labels")
    model = SCVIX(adata, encode_covariates=True)
    model.train(max_epochs=1)
    model.get_normalized_expression(n_samples=2)


def test_scvix_embedding():
    adata = synthetic_iid(batch_size=100)
    SCVIX.setup_anndata(adata, batch_key="batch", assay_key="batch", labels_key="labels")
    model = SCVIX(adata, batch_representation="embedding")
    model.train(max_epochs=1)
    model.get_normalized_expression(n_samples=2)


def test_scvix_layernorm():
    adata = synthetic_iid(batch_size=100)
    SCVIX.setup_anndata(adata, batch_key="batch", assay_key="batch", labels_key="labels")
    model = SCVIX(adata, conditional_norm=False, use_batch_norm="both", use_layer_norm="none")
    model.train(max_epochs=1)
    model.get_normalized_expression(n_samples=2)
    model = SCVIX(adata, conditional_norm=True, use_batch_norm="both", use_layer_norm="none")
    model.train(max_epochs=1)
    model.get_normalized_expression(n_samples=2)


def test_scvix_scarches_one_hot(save_path):
    # test transfer_anndata_setup + view
    adata1 = synthetic_iid()
    SCVIX.setup_anndata(adata1, batch_key="batch", assay_key="batch", labels_key="labels")
    model = SCVIX(adata1, batch_representation="one-hot")
    model.train(1, train_size=0.5)
    dir_path = os.path.join(save_path, "saved_model/")
    model.save(dir_path, overwrite=True)

    # adata2 has more genes and a perfect subset of adata1
    adata2 = synthetic_iid(n_genes=110)
    adata2.obs["batch"] = adata2.obs.batch.cat.rename_categories(["batch_2", "batch_3"])
    SCVIX.prepare_query_anndata(adata2, dir_path)
    SCVIX_query = SCVIX.load_query_data(adata2, dir_path)
    SCVIX_query.train(1, train_size=0.5, plan_kwargs={"weight_decay": 0.0})

    adata3 = SCVIX.prepare_query_anndata(adata2, dir_path, inplace=False)
    SCVIX_query2 = SCVIX.load_query_data(adata3, dir_path)
    SCVIX_query2.train(1, train_size=0.5, plan_kwargs={"weight_decay": 0.0})

    # adata4 has more genes and missing 10 genes from adata1
    adata4 = synthetic_iid(n_genes=110)
    new_var_names_init = [f"Random {i}" for i in range(10)]
    new_var_names = new_var_names_init + adata4.var_names[10:].to_list()
    adata4.var_names = new_var_names

    SCVIX.prepare_query_anndata(adata4, dir_path)
    # should be padded 0s
    assert np.sum(adata4[:, adata4.var_names[:10]].X) == 0
    np.testing.assert_equal(adata4.var_names[:10].to_numpy(), adata1.var_names[:10].to_numpy())
    SCVIX_query3 = SCVIX.load_query_data(adata4, dir_path)
    SCVIX_query3.train(1, train_size=0.5, plan_kwargs={"weight_decay": 0.0})

    adata5 = SCVIX.prepare_query_anndata(adata4, dir_path, inplace=False)
    SCVIX_query4 = SCVIX.load_query_data(adata5, dir_path)
    SCVIX_query4.train(1, train_size=0.5, plan_kwargs={"weight_decay": 0.0})


def test_scvix_scarches_embedding(save_path):
    # test transfer_anndata_setup + view
    adata1 = synthetic_iid()
    SCVIX.setup_anndata(adata1, batch_key="batch", assay_key="batch", labels_key="labels")
    model = SCVIX(adata1, batch_representation="embedding")
    model.train(1, train_size=0.5)
    dir_path = os.path.join(save_path, "saved_model/")
    model.save(dir_path, overwrite=True)

    # adata2 has more genes and a perfect subset of adata1
    adata2 = synthetic_iid(n_genes=110)
    adata2.obs["batch"] = adata2.obs.batch.cat.rename_categories(["batch_2", "batch_3"])
    SCVIX.prepare_query_anndata(adata2, dir_path)
    SCVIX_query = SCVIX.load_query_data(adata2, dir_path)
    SCVIX_query.train(1, train_size=0.5, plan_kwargs={"weight_decay": 0.0})

    adata3 = SCVIX.prepare_query_anndata(adata2, dir_path, inplace=False)
    SCVIX_query2 = SCVIX.load_query_data(adata3, dir_path)
    SCVIX_query2.train(1, train_size=0.5, plan_kwargs={"weight_decay": 0.0})

    # adata4 has more genes and missing 10 genes from adata1
    adata4 = synthetic_iid(n_genes=110)
    new_var_names_init = [f"Random {i}" for i in range(10)]
    new_var_names = new_var_names_init + adata4.var_names[10:].to_list()
    adata4.var_names = new_var_names

    SCVIX.prepare_query_anndata(adata4, dir_path)
    # should be padded 0s
    assert np.sum(adata4[:, adata4.var_names[:10]].X) == 0
    np.testing.assert_equal(adata4.var_names[:10].to_numpy(), adata1.var_names[:10].to_numpy())
    SCVIX_query3 = SCVIX.load_query_data(adata4, dir_path)
    SCVIX_query3.train(1, train_size=0.5, plan_kwargs={"weight_decay": 0.0})

    adata5 = SCVIX.prepare_query_anndata(adata4, dir_path, inplace=False)
    SCVIX_query4 = SCVIX.load_query_data(adata5, dir_path)
    SCVIX_query4.train(1, train_size=0.5, plan_kwargs={"weight_decay": 0.0})


def test_scvix_minified():
    n_samples = 200
    adata = synthetic_iid()
    SCVIX.setup_anndata(adata, batch_key="batch", assay_key="batch", labels_key="labels")
    model = SCVIX(adata, batch_representation="embedding", gene_likelihood="zinb")
    model.train(1, train_size=0.5)

    qzm, qzv = model.get_latent_representation(give_mean=False, return_dist=True)
    model.adata.obsm["X_latent_qzm"] = qzm
    model.adata.obsm["X_latent_qzv"] = qzv
    lib_size = np.squeeze(np.asarray(adata.X.sum(axis=-1)))

    scvi.settings.seed = 1
    params_orig = model.get_likelihood_parameters(n_samples=n_samples, give_mean=True)
    adata_orig = adata.copy()

    model.minify_adata()
    assert model.minified_data_type == ADATA_MINIFY_TYPE.LATENT_POSTERIOR
    assert model.adata_manager.registry is model.registry_

    assert not _is_minified(adata)
    assert adata is not model.adata

    orig_obs_df = adata_orig.obs
    orig_obs_df[BaseMinifiedModeModelClass._OBSERVED_LIB_SIZE_KEY] = lib_size
    assert model.adata.obs.equals(orig_obs_df)
    assert model.adata.var_names.equals(adata_orig.var_names)
    assert model.adata.var.equals(adata_orig.var)

    scvi.settings.seed = 1
    keys = ["mean", "dispersions", "dropout"]
    # The full-data path draws one unused encoder sample before the requested samples.
    # Draw and discard one cached sample so the stochastic latent and batch-embedding draws align.
    params_latent = model.get_likelihood_parameters(n_samples=n_samples + 1, give_mean=False)
    params_latent = {key: value[1:].mean(0) for key, value in params_latent.items()}
    for k in keys:
        assert params_latent[k].shape == params_orig[k].shape

    for k in keys:
        assert_approx_equal(params_latent[k], params_orig[k])
