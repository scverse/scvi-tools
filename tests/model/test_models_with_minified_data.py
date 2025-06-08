from __future__ import annotations

from typing import TYPE_CHECKING

import numpy as np
import pytest

import scvi
from scvi.data import synthetic_iid
from scvi.data._constants import _ADATA_MINIFY_TYPE_UNS_KEY, ADATA_MINIFY_TYPE
from scvi.data._utils import _is_minified
from scvi.model import SCANVI, SCVI
from scvi.model.base import BaseMinifiedModeModelClass

if TYPE_CHECKING:
    import numpy.typing as npt
    from anndata import AnnData


def prep_model(
    cls: BaseMinifiedModeModelClass = SCVI,
    layer: str | None = None,
    use_size_factor: bool = False,
    n_latent: int = 5,
) -> tuple[BaseMinifiedModeModelClass, AnnData, npt.NDArray, AnnData]:
    adata = synthetic_iid()
    counts = adata.X
    if use_size_factor:
        adata.obs["size_factor"] = np.random.randint(1, 5, size=(adata.shape[0],))
    if layer is not None:
        adata.layers[layer] = adata.X.copy()
        adata.X = np.zeros_like(adata.X)

    adata_before_setup = adata.copy()

    setup_kwargs = {
        "layer": layer,
        "batch_key": "batch",
        "labels_key": "labels",
        "size_factor_key": "size_factor" if use_size_factor else None,
    }
    if cls == SCANVI:
        setup_kwargs["unlabeled_category"] = "unknown"
    cls.setup_anndata(
        adata,
        **setup_kwargs,
    )

    model = cls(adata, n_latent=n_latent)
    model.train(1, check_val_every_n_epoch=1, train_size=0.5)

    lib_size = np.squeeze(np.asarray(counts.sum(axis=-1)))

    return model, adata, lib_size, adata_before_setup


def assert_approx_equal(a, b):
    # Allclose because on GPU, the values are not exactly the same
    # as some values are moved to cpu during data minification
    np.testing.assert_allclose(a, b, rtol=3e-1, atol=5e-1)


def run_test_for_model_with_minified_adata(
    cls: BaseMinifiedModeModelClass = SCVI,
    n_samples: int = 1,
    give_mean: bool = False,
    layer: str = None,
    use_size_factor: bool = False,
):
    model, adata, adata_lib_size, _ = prep_model(cls, layer, use_size_factor)

    scvi.settings.seed = 1
    qzm, qzv = model.get_latent_representation(give_mean=False, return_dist=True)
    model.adata.obsm["X_latent_qzm"] = qzm
    model.adata.obsm["X_latent_qzv"] = qzv
    scvi.settings.seed = 1
    params_orig = model.get_likelihood_parameters(n_samples=n_samples, give_mean=give_mean)
    adata_orig = adata.copy()

    model.minify_adata()
    assert model.minified_data_type == ADATA_MINIFY_TYPE.LATENT_POSTERIOR
    assert model.adata_manager.registry is model.registry_

    assert not _is_minified(adata)
    assert adata is not model.adata

    orig_obs_df = adata_orig.obs
    orig_obs_df[BaseMinifiedModeModelClass._OBSERVED_LIB_SIZE_KEY] = adata_lib_size
    assert model.adata.obs.equals(orig_obs_df)
    assert model.adata.var_names.equals(adata_orig.var_names)
    assert model.adata.var.equals(adata_orig.var)

    scvi.settings.seed = 1
    keys = ["mean", "dispersions", "dropout"]
    if n_samples == 1:
        params_latent = model.get_likelihood_parameters(n_samples=n_samples, give_mean=give_mean)
    else:
        # do this so that we generate the same sequence of random numbers in the
        # minified and non-minified cases (purely to get the tests to pass). this is
        # because in the non-minified case we sample once more (in the call to z_encoder
        # during inference)
        params_latent = model.get_likelihood_parameters(n_samples=n_samples, give_mean=give_mean)
    for k in keys:
        assert params_latent[k].shape == params_orig[k].shape

    for k in keys:
        assert_approx_equal(params_latent[k], params_orig[k])


def test_scvi_with_minified_adata_one_sample():
    run_test_for_model_with_minified_adata()
    run_test_for_model_with_minified_adata(layer="data_layer", use_size_factor=True)


def test_scvi_with_minified_adata_one_sample_with_layer():
    run_test_for_model_with_minified_adata(layer="data_layer")
    run_test_for_model_with_minified_adata(layer="data_layer", use_size_factor=True)


def test_scvi_with_minified_adata_n_samples():
    run_test_for_model_with_minified_adata(n_samples=400, give_mean=True)
    run_test_for_model_with_minified_adata(n_samples=400, give_mean=True, use_size_factor=True)


def test_scanvi_with_minified_adata_one_sample():
    run_test_for_model_with_minified_adata(SCANVI)
    run_test_for_model_with_minified_adata(SCANVI, use_size_factor=True)


def test_scanvi_with_minified_adata_one_sample_with_layer():
    run_test_for_model_with_minified_adata(SCANVI, layer="data_layer")
    run_test_for_model_with_minified_adata(SCANVI, layer="data_layer", use_size_factor=True)


def test_scanvi_with_minified_adata_n_samples():
    run_test_for_model_with_minified_adata(SCANVI, n_samples=400, give_mean=True)
    run_test_for_model_with_minified_adata(
        SCANVI, n_samples=400, give_mean=True, use_size_factor=True
    )


def test_scanvi_from_scvi(save_path):
    model, adata, _, adata_before_setup = prep_model(SCVI)
    qzm, qzv = model.get_latent_representation(give_mean=False, return_dist=True)
    adata.obsm["X_latent_qzm"] = qzm
    adata.obsm["X_latent_qzv"] = qzv
    model.minify_adata()
    assert model.minified_data_type == ADATA_MINIFY_TYPE.LATENT_POSTERIOR

    # for from_scvi_model to succeed:
    # 1. scvi_model.adata must not be minified, and
    # 2. adata (if different from scvi_model.adata) must not be minified

    with pytest.raises(ValueError) as e:
        scvi.model.SCANVI.from_scvi_model(model, "label_0")

    msg = (
        "We cannot use the given scVI model to initialize scANVI because it has minified adata. "
        "Keep counts when minifying model using minified_data_type="
        "'latent_posterior_parameters_with_counts'."
    )
    assert str(e.value) == msg

    # let's load scvi_model with a non-minified adata
    model.save(save_path, overwrite=True)
    loaded_model = SCVI.load(save_path, adata=adata_before_setup)

    adata2 = synthetic_iid()
    # just add this to pretend the data is minified
    adata2.uns[_ADATA_MINIFY_TYPE_UNS_KEY] = ADATA_MINIFY_TYPE.LATENT_POSTERIOR
    with pytest.raises(ValueError) as e:
        scvi.model.SCANVI.from_scvi_model(loaded_model, "label_0", adata=adata2)
    assert str(e.value) == "Please provide a non-minified `adata` to initialize scANVI."

    scanvi_model = scvi.model.SCANVI.from_scvi_model(loaded_model, "label_0")
    scanvi_model.train(1)


def test_scvi_with_minified_adata_get_normalized_expression():
    model, adata, _, _ = prep_model()

    scvi.settings.seed = 1
    qzm, qzv = model.get_latent_representation(give_mean=False, return_dist=True)
    model.adata.obsm["X_latent_qzm"] = qzm
    model.adata.obsm["X_latent_qzv"] = qzv

    scvi.settings.seed = 1
    exprs_orig = model.get_normalized_expression()

    model.minify_adata()
    assert model.minified_data_type == ADATA_MINIFY_TYPE.LATENT_POSTERIOR

    scvi.settings.seed = 1
    exprs_new = model.get_normalized_expression()
    assert exprs_new.shape == adata.shape

    np.testing.assert_allclose(exprs_new, exprs_orig, rtol=3e-1, atol=5e-1)


def test_scvi_with_minified_adata_get_normalized_expression_non_default_gene_list():
    model, adata, _, _ = prep_model()

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

    model.minify_adata()
    assert model.minified_data_type == ADATA_MINIFY_TYPE.LATENT_POSTERIOR

    scvi.settings.seed = 1
    # do this so that we generate the same sequence of random numbers in the
    # minified and non-minified cases (purely to get the tests to pass). this is
    # because in the non-minified case we sample once more (in the call to z_encoder
    # during inference)
    exprs_new = model.get_normalized_expression(
        gene_list=gl, n_samples=n_samples + 1, return_mean=False, library_size="latent"
    )
    exprs_new = exprs_new[1:].mean(0)

    assert exprs_new.shape == (adata.shape[0], 5)
    np.testing.assert_allclose(exprs_new, exprs_orig, rtol=3e-1, atol=5e-1)


def test_validate_unsupported_if_minified():
    model, _, _, _ = prep_model()

    qzm, qzv = model.get_latent_representation(give_mean=False, return_dist=True)
    model.adata.obsm["X_latent_qzm"] = qzm
    model.adata.obsm["X_latent_qzv"] = qzv

    model.minify_adata()
    assert model.minified_data_type == ADATA_MINIFY_TYPE.LATENT_POSTERIOR

    common_err_msg = "The {} function currently does not support minified data."

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

    with pytest.raises(ValueError) as e:
        model.train()
    assert str(e.value) == common_err_msg.format("VAE.loss")


def test_validate_supported_if_minified_keep_count():
    model, _, _, _ = prep_model()
    model2, _, _, _ = prep_model()

    qzm, qzv = model.get_latent_representation(give_mean=False, return_dist=True)
    model.adata.obsm["X_latent_qzm"] = qzm
    model.adata.obsm["X_latent_qzv"] = qzv

    model.minify_adata(minified_data_type="latent_posterior_parameters_with_counts")
    assert model.minified_data_type == ADATA_MINIFY_TYPE.LATENT_POSTERIOR_WITH_COUNTS
    assert model2.minified_data_type is None

    assert np.allclose(model2.get_elbo().cpu(), model.get_elbo().cpu(), rtol=5e-2)
    assert np.allclose(
        model2.get_reconstruction_error()["reconstruction_loss"].cpu(),
        model.get_reconstruction_error()["reconstruction_loss"].cpu(),
        rtol=5e-2,
    )
    assert np.allclose(model2.get_marginal_ll(), model.get_marginal_ll(), rtol=5e-2)

    model.train(1, check_val_every_n_epoch=1, train_size=0.5)
    model.train(
        1, check_val_every_n_epoch=1, train_size=0.5, plan_kwargs={"update_only_decoder": True}
    )
    scanvi_model = scvi.model.SCANVI.from_scvi_model(
        model, labels_key="labels", unlabeled_category="unknown"
    )
    scanvi_model.train()
    scanvi_model.train(
        1, check_val_every_n_epoch=1, train_size=0.5, plan_kwargs={"update_only_decoder": True}
    )


def test_scvi_with_minified_adata_save_then_load(save_path):
    # create a model and minify its adata, then save it and its adata.
    # Load it back up using the same (minified) adata. Validate that the
    # loaded model has the minified_data_type attribute set as expected.
    model, adata, _, _ = prep_model()

    scvi.settings.seed = 1
    qzm, qzv = model.get_latent_representation(give_mean=False, return_dist=True)
    model.adata.obsm["X_latent_qzm"] = qzm
    model.adata.obsm["X_latent_qzv"] = qzv

    scvi.settings.seed = 1
    params_orig = model.get_likelihood_parameters()

    model.minify_adata()
    assert model.minified_data_type == ADATA_MINIFY_TYPE.LATENT_POSTERIOR

    model.save(save_path, overwrite=True, save_anndata=True)
    # load saved model with saved (minified) adata
    loaded_model = SCVI.load(save_path)

    assert loaded_model.minified_data_type == ADATA_MINIFY_TYPE.LATENT_POSTERIOR

    scvi.settings.seed = 1
    params_latent = loaded_model.get_likelihood_parameters()
    assert params_latent["mean"].shape == adata.shape
    np.testing.assert_allclose(params_latent["mean"], params_orig["mean"], rtol=3e-1, atol=5e-1)


def test_scvi_with_minified_adata_save_then_load_with_non_minified_adata(save_path):
    # create a model and minify its adata, then save it and its adata.
    # Load it back up using a non-minified adata. Validate that the
    # loaded model does not has the minified_data_type attribute set.
    model, adata, _, adata_before_setup = prep_model()

    scvi.settings.seed = 1
    qzm, qzv = model.get_latent_representation(give_mean=False, return_dist=True)
    model.adata.obsm["X_latent_qzm"] = qzm
    model.adata.obsm["X_latent_qzv"] = qzv

    scvi.settings.seed = 1
    params_orig = model.get_likelihood_parameters()

    model.minify_adata()
    assert model.minified_data_type == ADATA_MINIFY_TYPE.LATENT_POSTERIOR

    model.save(save_path, overwrite=True, save_anndata=True)
    # load saved model with a non-minified adata
    loaded_model = SCVI.load(save_path, adata=adata_before_setup)

    assert loaded_model.minified_data_type is None

    scvi.settings.seed = 1
    params_new = loaded_model.get_likelihood_parameters()
    assert params_new["mean"].shape == adata.shape
    np.testing.assert_array_equal(params_new["mean"], params_orig["mean"])


def test_scvi_save_then_load_with_minified_adata(save_path):
    # create a model, then save it and its adata (non-minified).
    # Load it back up using a minified adata. Validate that this
    # fails, as expected because we don't have a way to validate
    # whether the minified-adata was set up correctly
    model, _, _, _ = prep_model()

    qzm, qzv = model.get_latent_representation(give_mean=False, return_dist=True)
    model.adata.obsm["X_latent_qzm"] = qzm
    model.adata.obsm["X_latent_qzv"] = qzv

    model.save(save_path, overwrite=True, save_anndata=True)

    model.minify_adata()
    assert model.minified_data_type == ADATA_MINIFY_TYPE.LATENT_POSTERIOR

    # loading this model with a minified adata is not allowed because
    # we don't have a way to validate whether the minified-adata was
    # set up correctly
    with pytest.raises(KeyError):
        SCVI.load(save_path, adata=model.adata)


def test_scvi_with_minified_adata_get_latent_representation():
    model, _, _, _ = prep_model()

    scvi.settings.seed = 1
    qzm, qzv = model.get_latent_representation(give_mean=False, return_dist=True)
    model.adata.obsm["X_latent_qzm"] = qzm
    model.adata.obsm["X_latent_qzv"] = qzv

    scvi.settings.seed = 1
    latent_repr_orig = model.get_latent_representation()

    model.minify_adata()
    assert model.minified_data_type == ADATA_MINIFY_TYPE.LATENT_POSTERIOR

    scvi.settings.seed = 1
    latent_repr_new = model.get_latent_representation()

    np.testing.assert_array_equal(latent_repr_new, latent_repr_orig)


def test_scvi_with_minified_adata_posterior_predictive_sample():
    model, _, _, _ = prep_model()

    scvi.settings.seed = 1
    qzm, qzv = model.get_latent_representation(give_mean=False, return_dist=True)
    model.adata.obsm["X_latent_qzm"] = qzm
    model.adata.obsm["X_latent_qzv"] = qzv

    scvi.settings.seed = 1
    sample_orig = model.posterior_predictive_sample(
        indices=[1, 2, 3], gene_list=["gene_1", "gene_2"]
    )

    model.minify_adata()
    assert model.minified_data_type == ADATA_MINIFY_TYPE.LATENT_POSTERIOR

    scvi.settings.seed = 1
    sample_new = model.posterior_predictive_sample(
        indices=[1, 2, 3], gene_list=["gene_1", "gene_2"]
    )
    assert sample_new.shape == (3, 2)

    np.testing.assert_array_equal(sample_new.todense(), sample_orig.todense())


def test_scvi_with_minified_adata_get_feature_correlation_matrix():
    model, _, _, _ = prep_model()

    scvi.settings.seed = 1
    qzm, qzv = model.get_latent_representation(give_mean=False, return_dist=True)
    model.adata.obsm["X_latent_qzm"] = qzm
    model.adata.obsm["X_latent_qzv"] = qzv

    scvi.settings.seed = 1
    fcm_orig = model.get_feature_correlation_matrix(
        correlation_type="pearson",
        n_samples=1,
        transform_batch=["batch_0", "batch_1"],
    )

    model.minify_adata()
    assert model.minified_data_type == ADATA_MINIFY_TYPE.LATENT_POSTERIOR

    scvi.settings.seed = 1
    fcm_new = model.get_feature_correlation_matrix(
        correlation_type="pearson",
        n_samples=1,
        transform_batch=["batch_0", "batch_1"],
    )

    assert_approx_equal(fcm_new, fcm_orig)
