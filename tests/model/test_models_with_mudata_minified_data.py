import numpy as np
import pytest

# from mudata import MuData
from scvi.data import synthetic_iid
from scvi.data._constants import ADATA_MINIFY_TYPE
from scvi.data._utils import _is_minified
from scvi.model import MULTIVI, TOTALVI

_TOTALVI_OBSERVED_LIB_SIZE = "_totalvi_observed_lib_size"
_MULTIVI_OBSERVED_LIB_SIZE = "_multivi_observed_lib_size"


def prep_model_mudata(cls=TOTALVI, use_size_factor=False):
    # create a synthetic dataset
    # adata = synthetic_iid()
    # protein_adata = synthetic_iid(n_genes=50)
    # mdata = MuData({"rna": adata, "protein_expression": protein_adata,
    #                 "accessibility": synthetic_iid()})
    mdata = synthetic_iid(return_mudata=True)
    if use_size_factor:
        mdata.obs["size_factor"] = np.random.randint(1, 5, size=(mdata.shape[0],))
    # if layer is not None:
    #     for mod in mdata.mod_names:
    #         mdata[mod].layers[layer] = mdata[mod].X.copy()
    #         mdata[mod].X = np.zeros_like(mdata[mod].X)
    mdata.var["n_counts"] = np.squeeze(
        np.concatenate(
            [
                np.asarray(np.sum(mdata["rna"].X, axis=0)),
                np.asarray(np.sum(mdata["protein_expression"].X, axis=0)),
                np.asarray(np.sum(mdata["accessibility"].X, axis=0)),
            ]
        )
    )
    mdata.varm["my_varm"] = np.random.negative_binomial(5, 0.3, size=(mdata.shape[1], 3))
    # mdata["rna"].layers["my_layer"] = np.ones_like(mdata["rna"].X)
    mdata_before_setup = mdata.copy()

    # run setup_anndata
    setup_kwargs = {
        "batch_key": "batch",
    }
    if use_size_factor:
        setup_kwargs["size_factor_key"] = "size_factor"

    if cls == TOTALVI:
        # create and train the model
        cls.setup_mudata(
            mdata,
            modalities={"rna_layer": "rna", "protein_layer": "protein_expression"},
            **setup_kwargs,
        )
        model = cls(mdata, n_latent=5)
    elif cls == MULTIVI:
        # create and train the model
        cls.setup_mudata(
            mdata,
            modalities={
                "rna_layer": "rna",
                "protein_layer": "protein_expression",
                "atac_layer": "accessibility",
            },
            **setup_kwargs,
        )
        model = cls(mdata, n_latent=5, n_genes=50, n_regions=50)
    else:
        raise ValueError("Bad Model name as input to test")
    model.train(1, check_val_every_n_epoch=1, train_size=0.5)

    # get the mdata lib size
    if cls == TOTALVI:
        mdata_lib_size = np.squeeze(np.asarray(mdata["rna"].X.sum(axis=1)))  # TOTALVI
    else:
        mdata_lib_size = np.squeeze(np.asarray(mdata["accessibility"].X.sum(axis=1)))  # MULTIVI
    assert (
        np.min(mdata_lib_size) > 0
    )  # make sure it's not all zeros and there are no negative values

    return model, mdata, mdata_lib_size, mdata_before_setup


def assert_approx_equal(a, b):
    # Allclose because on GPU, the values are not exactly the same
    # as some values are moved to cpu during data minification
    np.testing.assert_allclose(a, b, rtol=3e-1, atol=5e-1)


def run_test_for_model_with_minified_mudata(
    cls=TOTALVI,
    use_size_factor=False,
):
    model, mdata, mdata_lib_size, _ = prep_model_mudata(cls, use_size_factor)

    qzm, qzv = model.get_latent_representation(give_mean=False, return_dist=True)
    model.adata.obsm["X_latent_qzm"] = qzm
    model.adata.obsm["X_latent_qzv"] = qzv

    mdata_orig = mdata.copy()

    model.minify_mudata()
    assert model.minified_data_type == ADATA_MINIFY_TYPE.LATENT_POSTERIOR
    assert model.adata_manager.registry is model.registry_

    # make sure the original mdata we set up the model with was not changed
    assert mdata is not model.adata
    assert _is_minified(mdata) is False
    assert _is_minified(model.adata) is True

    assert mdata_orig["rna"].layers.keys() == model.adata["rna"].layers.keys()
    orig_obs_df = mdata_orig.obs
    obs_keys = _TOTALVI_OBSERVED_LIB_SIZE if cls == TOTALVI else _MULTIVI_OBSERVED_LIB_SIZE
    orig_obs_df[obs_keys] = mdata_lib_size
    assert model.adata.obs.equals(orig_obs_df)
    assert model.adata.var_names.equals(mdata_orig.var_names)
    assert model.adata.var.equals(mdata_orig.var)
    assert model.adata.varm.keys() == mdata_orig.varm.keys()
    np.testing.assert_array_equal(model.adata.varm["my_varm"], mdata_orig.varm["my_varm"])


@pytest.mark.parametrize("cls", [TOTALVI, MULTIVI])
@pytest.mark.parametrize("use_size_factor", [False, True])
def test_with_minified_mudata(cls, use_size_factor: bool):
    run_test_for_model_with_minified_mudata(cls=cls, use_size_factor=use_size_factor)


@pytest.mark.parametrize("cls", [TOTALVI, MULTIVI])
def test_with_minified_mdata_get_normalized_expression(cls):
    model, mdata, _, _ = prep_model_mudata(cls=cls)

    qzm, qzv = model.get_latent_representation(give_mean=False, return_dist=True)
    model.adata.obsm["X_latent_qzm"] = qzm
    model.adata.obsm["X_latent_qzv"] = qzv

    exprs_orig = model.get_normalized_expression()

    model.minify_mudata()
    assert model.minified_data_type == ADATA_MINIFY_TYPE.LATENT_POSTERIOR

    exprs_new = model.get_normalized_expression()
    if type(exprs_new) is tuple:
        for ii in range(len(exprs_new)):
            assert exprs_new[ii].shape == mdata[mdata.mod_names[ii]].shape
        for ii in range(len(exprs_new)):
            np.testing.assert_array_equal(exprs_new[ii], exprs_orig[ii])
    else:
        assert exprs_new.shape == mdata[mdata.mod_names].shape
        np.testing.assert_array_equal(exprs_new, exprs_orig)


@pytest.mark.parametrize("cls", [TOTALVI, MULTIVI])
def test_with_minified_mdata_get_normalized_expression_non_default_gene_list(cls):
    model, mdata, _, _ = prep_model_mudata(cls=cls)

    # non-default gene list and n_samples > 1
    gl = mdata.var_names[:5].to_list()
    n_samples = 10

    qzm, qzv = model.get_latent_representation(give_mean=False, return_dist=True)
    model.adata.obsm["X_latent_qzm"] = qzm
    model.adata.obsm["X_latent_qzv"] = qzv

    exprs_orig = model.get_normalized_expression(
        gene_list=gl, n_samples=n_samples, library_size="latent"
    )

    model.minify_mudata()
    assert model.minified_data_type == ADATA_MINIFY_TYPE.LATENT_POSTERIOR

    # do this so that we generate the same sequence of random numbers in the
    # minified and non-minified cases (purely to get the tests to pass). this is
    # because in the non-minified case we sample once more (in the call to z_encoder
    # during inference)
    exprs_new = model.get_normalized_expression(
        gene_list=gl, n_samples=n_samples + 1, return_mean=False, library_size="latent"
    )
    exprs_new = exprs_new[0][:, :, 1:].mean(2)

    assert exprs_new.shape == (mdata.shape[0], 5)
    np.testing.assert_allclose(exprs_new, exprs_orig[0], rtol=3e-1, atol=5e-1)


@pytest.mark.parametrize("cls", [TOTALVI, MULTIVI])
def test_validate_unsupported_if_minified(cls):
    model, _, _, _ = prep_model_mudata(cls=cls)

    qzm, qzv = model.get_latent_representation(give_mean=False, return_dist=True)
    model.adata.obsm["X_latent_qzm"] = qzm
    model.adata.obsm["X_latent_qzv"] = qzv

    model.minify_mudata()
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

    if cls != TOTALVI:
        with pytest.raises(AttributeError) as e:
            model.get_latent_library_size()
        assert str(e.value) == "'MULTIVI' object has no attribute 'get_latent_library_size'"


@pytest.mark.parametrize("cls", [TOTALVI, MULTIVI])
def test_with_minified_mdata_save_then_load(cls, save_path):
    # create a model and minify its mdata, then save it and its mdata.
    # Load it back up using the same (minified) mdata. Validate that the
    # loaded model has the minified_data_type attribute set as expected.
    model, mdata, _, _ = prep_model_mudata(cls=cls)

    qzm, qzv = model.get_latent_representation(give_mean=False, return_dist=True)
    model.adata.obsm["X_latent_qzm"] = qzm
    model.adata.obsm["X_latent_qzv"] = qzv

    model.minify_mudata()
    assert model.minified_data_type == ADATA_MINIFY_TYPE.LATENT_POSTERIOR

    model.save(save_path, overwrite=True, save_anndata=True, legacy_mudata_format=True)
    model.view_setup_args(save_path)
    # load saved model with saved (minified) mdata
    loaded_model = cls.load(save_path, adata=mdata)

    assert loaded_model.minified_data_type == ADATA_MINIFY_TYPE.LATENT_POSTERIOR


@pytest.mark.parametrize("cls", [TOTALVI, MULTIVI])
def test_with_minified_mdata_save_then_load_with_non_minified_mdata(cls, save_path):
    # create a model and minify its mdata, then save it and its mdata.
    # Load it back up using a non-minified mdata. Validate that the
    # loaded model does not has the minified_data_type attribute set.
    model, mdata, _, mdata_before_setup = prep_model_mudata(cls=cls)

    qzm, qzv = model.get_latent_representation(give_mean=False, return_dist=True)
    model.adata.obsm["X_latent_qzm"] = qzm
    model.adata.obsm["X_latent_qzv"] = qzv

    model.minify_mudata()
    assert model.minified_data_type == ADATA_MINIFY_TYPE.LATENT_POSTERIOR

    model.save(save_path, overwrite=True, save_anndata=False, legacy_mudata_format=True)
    # load saved model with a non-minified mdata
    loaded_model = cls.load(save_path, adata=mdata_before_setup)

    assert loaded_model.minified_data_type is None


@pytest.mark.parametrize("cls", [TOTALVI, MULTIVI])
def test_save_then_load_with_minified_mdata(cls, save_path):
    # create a model, then save it and its mdata (non-minified).
    # Load it back up using a minified mdata. Validate that this
    # fails, as expected because we don't have a way to validate
    # whether the minified-mdata was set up correctly
    model, _, _, _ = prep_model_mudata(cls=cls)

    qzm, qzv = model.get_latent_representation(give_mean=False, return_dist=True)
    model.adata.obsm["X_latent_qzm"] = qzm
    model.adata.obsm["X_latent_qzv"] = qzv

    model.save(save_path, overwrite=True, save_anndata=False, legacy_mudata_format=True)

    model.minify_mudata()
    assert model.minified_data_type == ADATA_MINIFY_TYPE.LATENT_POSTERIOR

    # loading this model with a minified mdata is not allowed because
    # we don't have a way to validate whether the minified-mdata was
    # set up correctly
    with pytest.raises(KeyError):
        cls.load(save_path, adata=model.adata)


@pytest.mark.parametrize("cls", [TOTALVI, MULTIVI])
def test_with_minified_mdata_get_latent_representation(cls):
    model, _, _, _ = prep_model_mudata(cls=cls)

    qzm, qzv = model.get_latent_representation(give_mean=False, return_dist=True)
    model.adata.obsm["X_latent_qzm"] = qzm
    model.adata.obsm["X_latent_qzv"] = qzv

    latent_repr_orig = model.get_latent_representation()

    model.minify_mudata()
    assert model.minified_data_type == ADATA_MINIFY_TYPE.LATENT_POSTERIOR

    latent_repr_new = model.get_latent_representation()

    np.testing.assert_array_equal(latent_repr_new, latent_repr_orig)


@pytest.mark.parametrize("cls", [TOTALVI])
def test_with_minified_mdata_posterior_predictive_sample(cls):
    model, _, _, _ = prep_model_mudata(cls=cls)

    qzm, qzv = model.get_latent_representation(give_mean=False, return_dist=True)
    model.adata.obsm["X_latent_qzm"] = qzm
    model.adata.obsm["X_latent_qzv"] = qzv

    sample_orig = model.posterior_predictive_sample()

    model.minify_mudata()
    assert model.minified_data_type == ADATA_MINIFY_TYPE.LATENT_POSTERIOR

    sample_new = model.posterior_predictive_sample()
    # assert sample_new.shape == (3, 2)

    np.testing.assert_array_equal(sample_new, sample_orig)


@pytest.mark.parametrize("cls", [TOTALVI])
def test_with_minified_mdata_get_feature_correlation_matrix(cls):
    model, _, _, _ = prep_model_mudata(cls=cls)

    qzm, qzv = model.get_latent_representation(give_mean=False, return_dist=True)
    model.adata.obsm["X_latent_qzm"] = qzm
    model.adata.obsm["X_latent_qzv"] = qzv

    fcm_orig = model.get_feature_correlation_matrix(
        correlation_type="spearman",
        transform_batch=["batch_0", "batch_1"],
    )

    model.minify_mudata()
    assert model.minified_data_type == ADATA_MINIFY_TYPE.LATENT_POSTERIOR

    fcm_new = model.get_feature_correlation_matrix(
        correlation_type="spearman",
        transform_batch=["batch_0", "batch_1"],
    )

    assert_approx_equal(fcm_new, fcm_orig)
