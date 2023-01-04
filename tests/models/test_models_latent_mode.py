from copy import deepcopy

import numpy as np
import pytest

import scvi
from scvi.data import synthetic_iid
from scvi.model import SCVI


def prep_model(layer=None):
    adata = synthetic_iid()
    adata.obs["size_factor"] = np.random.randint(1, 5, size=(adata.shape[0],))
    if layer is not None:
        adata.layers[layer] = adata.X.copy()
        adata.X = np.zeros_like(adata.X)
    adata.var["n_counts"] = np.squeeze(np.asarray(np.sum(adata.X, axis=0)))
    adata.varm["my_varm"] = np.random.negative_binomial(
        5, 0.3, size=(adata.shape[1], 3)
    )
    adata.layers["my_layer"] = np.ones_like(adata.X)
    SCVI.setup_anndata(
        adata,
        layer=layer,
        batch_key="batch",
        labels_key="labels",
        size_factor_key="size_factor",
    )
    model = SCVI(adata, n_latent=5)
    model.train(1, check_val_every_n_epoch=1, train_size=0.5)

    return model, adata


def run_test_scvi_latent_mode_dist(
    n_samples: int = 1, give_mean: bool = False, layer: str = None
):
    model, adata = prep_model(layer)

    scvi.settings.seed = 1
    qzm, qzv = model.get_latent_representation(give_mean=False, return_dist=True)
    adata.obsm["X_latent_qzm"] = qzm
    adata.obsm["X_latent_qzv"] = qzv

    scvi.settings.seed = 1
    params_orig = model.get_likelihood_parameters(
        n_samples=n_samples, give_mean=give_mean
    )
    model_orig = deepcopy(model)

    model.to_latent_mode()

    assert model.latent_data_type == "posterior_parameters"

    assert model_orig.adata.layers.keys() == model.adata.layers.keys()
    assert model.adata.obs.equals(model_orig.adata.obs)
    assert model.adata.var_names.equals(model_orig.adata.var_names)
    assert model.adata.var.equals(model_orig.adata.var)
    assert model.adata.varm.keys() == model_orig.adata.varm.keys()
    np.testing.assert_array_equal(
        model.adata.varm["my_varm"], model_orig.adata.varm["my_varm"]
    )

    scvi.settings.seed = 1
    keys = ["mean", "dispersions", "dropout"]
    if n_samples == 1:
        params_latent = model.get_likelihood_parameters(
            n_samples=n_samples, give_mean=give_mean
        )
    else:
        # do this so that we generate the same sequence of random numbers in the
        # latent and non latent cases (purely to get the tests to pass). this is
        # because in the non latent case we sample once more (in the call to z_encoder
        # during inference)
        params_latent = model.get_likelihood_parameters(
            n_samples=n_samples + 1, give_mean=False
        )
        for k in keys:
            params_latent[k] = params_latent[k][1:].mean(0)
    for k in keys:
        assert params_latent[k].shape == adata.shape

    for k in keys:
        # Allclose because on GPU, the values are not exactly the same
        # as latents are moved to cpu in latent mode
        np.testing.assert_allclose(
            params_latent[k], params_orig[k], rtol=3e-1, atol=5e-1
        )


def test_scvi_latent_mode_dist_one_sample():
    run_test_scvi_latent_mode_dist()


def test_scvi_latent_mode_dist_one_sample_with_layer():
    run_test_scvi_latent_mode_dist(layer="data_layer")


def test_scvi_latent_mode_dist_n_samples():
    run_test_scvi_latent_mode_dist(n_samples=10, give_mean=True)


def test_scvi_latent_mode_get_normalized_expression():
    model, adata = prep_model()

    scvi.settings.seed = 1
    qzm, qzv = model.get_latent_representation(give_mean=False, return_dist=True)
    adata.obsm["X_latent_qzm"] = qzm
    adata.obsm["X_latent_qzv"] = qzv

    scvi.settings.seed = 1
    exprs_orig = model.get_normalized_expression()

    model.to_latent_mode()

    scvi.settings.seed = 1
    exprs_latent = model.get_normalized_expression()
    assert exprs_latent.shape == adata.shape

    np.testing.assert_array_equal(exprs_latent, exprs_orig)


def test_scvi_latent_mode_get_normalized_expression_non_default_gene_list():
    model, adata = prep_model()

    # non-default gene list and n_samples > 1
    gl = adata.var_names[:5].to_list()
    n_samples = 10

    scvi.settings.seed = 1
    qzm, qzv = model.get_latent_representation(give_mean=False, return_dist=True)
    model.adata.obsm["X_latent_qzm"] = qzm
    model.adata.obsm["X_latent_qzv"] = qzv

    scvi.settings.seed = 1
    exprs_orig = model.get_normalized_expression(
        gene_list=gl, n_samples=n_samples, library_size="latent"
    )

    model.to_latent_mode()

    scvi.settings.seed = 1
    # do this so that we generate the same sequence of random numbers in the
    # latent and non latent cases (purely to get the tests to pass). this is
    # because in the non latent case we sample once more (in the call to z_encoder
    # during inference)
    exprs_latent = model.get_normalized_expression(
        gene_list=gl, n_samples=n_samples + 1, return_mean=False, library_size="latent"
    )
    exprs_latent = exprs_latent[1:].mean(0)

    assert exprs_latent.shape == (adata.shape[0], 5)
    np.testing.assert_allclose(exprs_latent, exprs_orig, rtol=3e-1, atol=5e-1)


def test_latent_mode_validate_unsupported():
    model, adata = prep_model()

    qzm, qzv = model.get_latent_representation(give_mean=False, return_dist=True)
    adata.obsm["X_latent_qzm"] = qzm
    adata.obsm["X_latent_qzv"] = qzv

    model.to_latent_mode()

    common_err_msg = "Latent mode currently not supported for the {} function."

    with pytest.raises(ValueError) as e:
        model.get_elbo()
    assert str(e.value) == common_err_msg.format("VAEMixin.get_elbo")

    with pytest.raises(ValueError) as e:
        model.get_reconstruction_error()
    assert str(e.value) == common_err_msg.format("VAEMixin.get_reconstruction_error")

    with pytest.raises(ValueError) as e:
        model.get_marginal_ll()
    assert str(e.value) == common_err_msg.format("VAEMixin.get_marginal_ll")

    with pytest.raises(ValueError) as e:
        model.get_latent_library_size()
    assert str(e.value) == common_err_msg.format("RNASeqMixin.get_latent_library_size")


def test_scvi_latent_mode_save_load_latent(save_path):
    model, adata = prep_model()

    scvi.settings.seed = 1
    qzm, qzv = model.get_latent_representation(give_mean=False, return_dist=True)
    adata.obsm["X_latent_qzm"] = qzm
    adata.obsm["X_latent_qzv"] = qzv

    scvi.settings.seed = 1
    params_orig = model.get_likelihood_parameters()

    model.to_latent_mode()

    model.save(save_path, overwrite=True, save_anndata=True)
    # load saved latent model with saved latent adata
    loaded_model = SCVI.load(save_path)

    scvi.settings.seed = 1
    params_latent = loaded_model.get_likelihood_parameters()
    assert params_latent["mean"].shape == adata.shape
    np.testing.assert_array_equal(params_latent["mean"], params_orig["mean"])


def test_scvi_latent_mode_save_load_latent_to_non_latent(save_path):
    model, adata = prep_model()

    scvi.settings.seed = 1
    qzm, qzv = model.get_latent_representation(give_mean=False, return_dist=True)
    adata.obsm["X_latent_qzm"] = qzm
    adata.obsm["X_latent_qzv"] = qzv

    scvi.settings.seed = 1
    params_orig = model.get_likelihood_parameters()
    model_orig = deepcopy(model)

    model.to_latent_mode()

    model.save(save_path, overwrite=True, save_anndata=True)
    # load saved latent model with non-latent adata
    loaded_model = SCVI.load(save_path, adata=model_orig.adata)

    scvi.settings.seed = 1
    params_new = loaded_model.get_likelihood_parameters()
    assert params_new["mean"].shape == adata.shape
    np.testing.assert_array_equal(params_new["mean"], params_orig["mean"])


def test_scvi_latent_mode_save_load_non_latent_to_latent(save_path):
    model, adata = prep_model()

    qzm, qzv = model.get_latent_representation(give_mean=False, return_dist=True)
    adata.obsm["X_latent_qzm"] = qzm
    adata.obsm["X_latent_qzv"] = qzv

    model_orig = deepcopy(model)
    model.to_latent_mode()

    model_orig.save(save_path, overwrite=True, save_anndata=True)

    # loading saved non-latent model with latent adata is not allowed
    # because we don't have a way to validate the correctness of the
    # latent data setup
    with pytest.raises(KeyError):
        SCVI.load(save_path, adata=model.adata)


def test_scvi_latent_mode_get_latent_representation():
    model, adata = prep_model()

    scvi.settings.seed = 1
    qzm, qzv = model.get_latent_representation(give_mean=False, return_dist=True)
    adata.obsm["X_latent_qzm"] = qzm
    adata.obsm["X_latent_qzv"] = qzv

    scvi.settings.seed = 1
    latent_repr_orig = model.get_latent_representation()

    model.to_latent_mode()

    scvi.settings.seed = 1
    latent_repr_latent = model.get_latent_representation()

    np.testing.assert_array_equal(latent_repr_latent, latent_repr_orig)


def test_scvi_latent_mode_posterior_predictive_sample():
    model, adata = prep_model()

    scvi.settings.seed = 1
    qzm, qzv = model.get_latent_representation(give_mean=False, return_dist=True)
    adata.obsm["X_latent_qzm"] = qzm
    adata.obsm["X_latent_qzv"] = qzv

    scvi.settings.seed = 1
    sample_orig = model.posterior_predictive_sample(
        indices=[1, 2, 3], gene_list=["1", "2"]
    )

    model.to_latent_mode()

    scvi.settings.seed = 1
    sample_latent = model.posterior_predictive_sample(
        indices=[1, 2, 3], gene_list=["1", "2"]
    )
    assert sample_latent.shape == (3, 2)

    np.testing.assert_array_equal(sample_latent, sample_orig)


def test_scvi_latent_mode_get_feature_correlation_matrix():
    model, adata = prep_model()

    scvi.settings.seed = 1
    qzm, qzv = model.get_latent_representation(give_mean=False, return_dist=True)
    adata.obsm["X_latent_qzm"] = qzm
    adata.obsm["X_latent_qzv"] = qzv

    scvi.settings.seed = 1
    fcm_orig = model.get_feature_correlation_matrix(
        correlation_type="pearson",
        n_samples=1,
        transform_batch=["batch_0", "batch_1"],
    )

    model.to_latent_mode()

    scvi.settings.seed = 1
    fcm_latent = model.get_feature_correlation_matrix(
        correlation_type="pearson",
        n_samples=1,
        transform_batch=["batch_0", "batch_1"],
    )

    np.testing.assert_array_equal(fcm_latent, fcm_orig)
