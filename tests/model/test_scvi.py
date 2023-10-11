import inspect
import os
import pickle
import tarfile
from unittest import mock

import anndata
import numpy as np
import pytest
import torch
from lightning.pytorch.callbacks import LearningRateMonitor
from scipy.sparse import csr_matrix
from torch.nn import Softplus

import scvi
from scvi.data import _constants, synthetic_iid
from scvi.data._compat import LEGACY_REGISTRY_KEY_MAP, registry_from_setup_dict
from scvi.data._download import _download
from scvi.model import SCVI
from scvi.model.utils import mde
from scvi.utils import attrdict

LEGACY_REGISTRY_KEYS = set(LEGACY_REGISTRY_KEY_MAP.values())
LEGACY_SETUP_DICT = {
    "scvi_version": "0.0.0",
    "categorical_mappings": {
        "_scvi_batch": {
            "original_key": "testbatch",
            "mapping": np.array(["batch_0", "batch_1"], dtype=object),
        },
        "_scvi_labels": {
            "original_key": "testlabels",
            "mapping": np.array(["label_0", "label_1", "label_2"], dtype=object),
        },
    },
    "extra_categoricals": {
        "mappings": {
            "cat1": np.array([0, 1, 2, 3, 4]),
            "cat2": np.array([0, 1, 2, 3, 4]),
        },
        "keys": ["cat1", "cat2"],
        "n_cats_per_key": [5, 5],
    },
    "extra_continuous_keys": np.array(["cont1", "cont2"], dtype=object),
    "data_registry": {
        "X": {"attr_name": "X", "attr_key": None},
        "batch_indices": {"attr_name": "obs", "attr_key": "_scvi_batch"},
        "labels": {"attr_name": "obs", "attr_key": "_scvi_labels"},
        "cat_covs": {
            "attr_name": "obsm",
            "attr_key": "_scvi_extra_categoricals",
        },
        "cont_covs": {
            "attr_name": "obsm",
            "attr_key": "_scvi_extra_continuous",
        },
    },
    "summary_stats": {
        "n_batch": 2,
        "n_cells": 400,
        "n_vars": 100,
        "n_labels": 3,
        "n_proteins": 0,
        "n_continuous_covs": 2,
    },
}


def test_saving_and_loading(save_path):
    def legacy_save(
        model,
        dir_path,
        prefix=None,
        overwrite=False,
        save_anndata=False,
        **anndata_write_kwargs,
    ):
        if not os.path.exists(dir_path) or overwrite:
            os.makedirs(dir_path, exist_ok=overwrite)
        else:
            raise ValueError(
                f"{dir_path} already exists. Please provide an unexisting directory for saving."
            )

        file_name_prefix = prefix or ""

        if save_anndata:
            model.adata.write(
                os.path.join(dir_path, f"{file_name_prefix}adata.h5ad"),
                **anndata_write_kwargs,
            )

        model_save_path = os.path.join(dir_path, f"{file_name_prefix}model_params.pt")
        attr_save_path = os.path.join(dir_path, f"{file_name_prefix}attr.pkl")
        varnames_save_path = os.path.join(dir_path, f"{file_name_prefix}var_names.csv")

        torch.save(model.module.state_dict(), model_save_path)

        var_names = model.adata.var_names.astype(str)
        var_names = var_names.to_numpy()
        np.savetxt(varnames_save_path, var_names, fmt="%s")

        # get all the user attributes
        user_attributes = model._get_user_attributes()
        # only save the public attributes with _ at the very end
        user_attributes = {a[0]: a[1] for a in user_attributes if a[0][-1] == "_"}

        with open(attr_save_path, "wb") as f:
            pickle.dump(user_attributes, f)

    def test_save_load_model(cls, adata, save_path, prefix=None):
        cls.setup_anndata(adata, batch_key="batch", labels_key="labels")
        model = cls(adata, latent_distribution="normal")
        model.train(1, train_size=0.2)
        z1 = model.get_latent_representation(adata)
        test_idx1 = model.validation_indices
        model.save(save_path, overwrite=True, save_anndata=True, prefix=prefix)
        model.view_setup_args(save_path, prefix=prefix)
        model = cls.load(save_path, prefix=prefix)
        model.get_latent_representation()

        # Load with mismatched genes.
        tmp_adata = synthetic_iid(
            n_genes=200,
        )
        with pytest.raises(ValueError):
            cls.load(save_path, adata=tmp_adata, prefix=prefix)

        # Load with different batches.
        tmp_adata = synthetic_iid()
        tmp_adata.obs["batch"] = tmp_adata.obs["batch"].cat.rename_categories(
            ["batch_2", "batch_3"]
        )
        with pytest.raises(ValueError):
            cls.load(save_path, adata=tmp_adata, prefix=prefix)

        model = cls.load(save_path, adata=adata, prefix=prefix)
        assert "batch" in model.adata_manager.data_registry
        assert model.adata_manager.data_registry.batch == attrdict(
            {"attr_name": "obs", "attr_key": "_scvi_batch"}
        )

        z2 = model.get_latent_representation()
        test_idx2 = model.validation_indices
        np.testing.assert_array_equal(z1, z2)
        np.testing.assert_array_equal(test_idx1, test_idx2)
        assert model.is_trained is True

        # Test legacy loading
        legacy_save_path = os.path.join(save_path, "legacy/")
        legacy_save(
            model, legacy_save_path, overwrite=True, save_anndata=True, prefix=prefix
        )
        with pytest.raises(ValueError):
            cls.load(legacy_save_path, adata=adata, prefix=prefix)
        cls.convert_legacy_save(
            legacy_save_path,
            legacy_save_path,
            overwrite=True,
            prefix=prefix,
        )
        m = cls.load(legacy_save_path, adata=adata, prefix=prefix)
        m.train(1)

    save_path = os.path.join(save_path, "tmp")
    adata = synthetic_iid()

    test_save_load_model(SCVI, adata, save_path, prefix=f"{SCVI.__name__}_")


def test_scvi(n_latent: int = 5):
    adata = synthetic_iid()
    adata.obs["size_factor"] = np.random.randint(1, 5, size=(adata.shape[0],))
    SCVI.setup_anndata(
        adata,
        batch_key="batch",
        labels_key="labels",
        size_factor_key="size_factor",
    )
    model = SCVI(adata, n_latent=n_latent)
    model.train(1, check_val_every_n_epoch=1, train_size=0.5)

    # test mde
    mde(model.get_latent_representation())

    # Test with observed lib size.
    adata = synthetic_iid()
    SCVI.setup_anndata(
        adata,
        batch_key="batch",
        labels_key="labels",
    )
    model = SCVI(adata, n_latent=n_latent)
    model.train(1, check_val_every_n_epoch=1, train_size=0.5)

    # Test without observed lib size.
    model = SCVI(
        adata, n_latent=n_latent, var_activation=Softplus(), use_observed_lib_size=False
    )
    model.train(1, check_val_every_n_epoch=1, train_size=0.5)
    model.train(1, check_val_every_n_epoch=1, train_size=0.5)

    # tests __repr__
    print(model)
    # test view_anndata_setup
    model.view_anndata_setup()
    model.view_anndata_setup(hide_state_registries=True)

    assert model.is_trained is True
    z = model.get_latent_representation()
    assert z.shape == (adata.shape[0], n_latent)
    assert len(model.history["elbo_train"]) == 2
    model.get_elbo()
    model.get_marginal_ll(n_mc_samples=3)
    model.get_reconstruction_error()
    model.get_normalized_expression(transform_batch="batch_1")
    model.get_normalized_expression(n_samples=2)

    adata2 = synthetic_iid()
    # test view_anndata_setup with different anndata before transfer setup
    with pytest.raises(ValueError):
        model.view_anndata_setup(adata=adata2)
        model.view_anndata_setup(adata=adata2, hide_state_registries=True)
    # test get methods with different anndata
    model.get_elbo(adata2)
    model.get_marginal_ll(adata2, n_mc_samples=3)
    model.get_reconstruction_error(adata2)
    latent = model.get_latent_representation(adata2, indices=[1, 2, 3])
    assert latent.shape == (3, n_latent)
    denoised = model.get_normalized_expression(adata2)
    assert denoised.shape == adata.shape
    # test view_anndata_setup with different anndata after transfer setup
    model.view_anndata_setup(adata=adata2)
    model.view_anndata_setup(adata=adata2, hide_state_registries=True)

    denoised = model.get_normalized_expression(
        adata2, indices=[1, 2, 3], transform_batch="batch_1"
    )
    denoised = model.get_normalized_expression(
        adata2, indices=[1, 2, 3], transform_batch=["batch_0", "batch_1"]
    )
    assert denoised.shape == (3, adata2.n_vars)
    sample = model.posterior_predictive_sample(adata2)
    assert sample.shape == adata2.shape
    sample = model.posterior_predictive_sample(
        adata2, indices=[1, 2, 3], gene_list=["1", "2"]
    )
    assert sample.shape == (3, 2)
    sample = model.posterior_predictive_sample(
        adata2, indices=[1, 2, 3], gene_list=["1", "2"], n_samples=3
    )
    assert sample.shape == (3, 2, 3)

    model.get_feature_correlation_matrix(correlation_type="pearson")
    model.get_feature_correlation_matrix(
        adata2,
        indices=[1, 2, 3],
        correlation_type="spearman",
        rna_size_factor=500,
        n_samples=5,
    )
    model.get_feature_correlation_matrix(
        adata2,
        indices=[1, 2, 3],
        correlation_type="spearman",
        rna_size_factor=500,
        n_samples=5,
        transform_batch=["batch_0", "batch_1"],
    )
    params = model.get_likelihood_parameters()
    assert params["mean"].shape == adata.shape
    assert (
        params["mean"].shape == params["dispersions"].shape == params["dropout"].shape
    )
    params = model.get_likelihood_parameters(adata2, indices=[1, 2, 3])
    assert params["mean"].shape == (3, adata.n_vars)
    params = model.get_likelihood_parameters(
        adata2, indices=[1, 2, 3], n_samples=3, give_mean=True
    )
    assert params["mean"].shape == (3, adata.n_vars)
    model.get_latent_library_size()
    model.get_latent_library_size(adata2, indices=[1, 2, 3])

    # test transfer_anndata_setup
    adata2 = synthetic_iid()
    model._validate_anndata(adata2)
    model.get_elbo(adata2)

    # test automatic transfer_anndata_setup on a view
    adata = synthetic_iid()
    SCVI.setup_anndata(
        adata,
        batch_key="batch",
        labels_key="labels",
    )
    model = SCVI(adata)
    adata2 = synthetic_iid()
    model.get_elbo(adata2[:10])

    # test automatic transfer_anndata_setup on a copy
    adata = synthetic_iid()
    SCVI.setup_anndata(
        adata,
        batch_key="batch",
        labels_key="labels",
    )
    model = SCVI(adata)
    model.train(1, train_size=0.5)
    adata2 = adata.copy()
    model.get_elbo(adata2)
    assert adata.uns[_constants._SCVI_UUID_KEY] != adata2.uns[_constants._SCVI_UUID_KEY]

    # test mismatched categories raises ValueError
    adata2 = synthetic_iid()
    adata2.obs.labels = adata2.obs.labels.cat.rename_categories(["a", "b", "c"])
    with pytest.raises(ValueError):
        model.get_elbo(adata2)

    # test differential expression
    model.differential_expression(groupby="labels", group1="label_1")
    model.differential_expression(
        groupby="labels", group1="label_1", weights="importance"
    )
    model.differential_expression(
        groupby="labels", group1="label_1", group2="label_2", mode="change"
    )
    model.differential_expression(groupby="labels")
    model.differential_expression(idx1=[0, 1, 2], idx2=[3, 4, 5])
    model.differential_expression(idx1=[0, 1, 2], idx2=[3, 4, 5], weights="importance")
    model.differential_expression(idx1=[0, 1, 2], idx2=[3, 4, 5])
    model.differential_expression(idx1=[0, 1, 2])

    model2 = SCVI(adata, use_observed_lib_size=False)
    model2.train(1)
    model2.differential_expression(idx1=[0, 1, 2], idx2=[3, 4, 5])
    model2.differential_expression(idx1=[0, 1, 2], idx2=[3, 4, 5], weights="importance")

    # transform batch works with all different types
    a = synthetic_iid()
    batch = np.zeros(a.n_obs)
    batch[:64] += 1
    a.obs["batch"] = batch
    SCVI.setup_anndata(
        a,
        batch_key="batch",
    )
    m = SCVI(a)
    m.train(1, train_size=0.5)
    m.get_normalized_expression(transform_batch=1)
    m.get_normalized_expression(transform_batch=[0, 1])

    # test get_likelihood_parameters() when dispersion=='gene-cell'
    model = SCVI(adata, dispersion="gene-cell")
    model.get_likelihood_parameters()
    model.get_likelihood_parameters(indices=np.arange(10))
    model.get_likelihood_parameters(n_samples=10)
    model.get_likelihood_parameters(n_samples=10, indices=np.arange(10))

    # test get_likelihood_parameters() when gene_likelihood!='zinb'
    model = SCVI(adata, gene_likelihood="nb")
    model.get_likelihood_parameters()

    # test different gene_likelihoods
    for gene_likelihood in ["zinb", "nb", "poisson"]:
        model = SCVI(adata, gene_likelihood=gene_likelihood)
        model.train(1, check_val_every_n_epoch=1, train_size=0.5)
        model.posterior_predictive_sample()
        model.get_latent_representation()
        model.get_normalized_expression()

    # test train callbacks work
    a = synthetic_iid()
    SCVI.setup_anndata(
        a,
        batch_key="batch",
        labels_key="labels",
    )
    m = SCVI(a)
    lr_monitor = LearningRateMonitor()
    m.train(
        callbacks=[lr_monitor],
        max_epochs=10,
        check_val_every_n_epoch=1,
        log_every_n_steps=1,
        plan_kwargs={"reduce_lr_on_plateau": True},
    )
    assert "lr-Adam" in m.history.keys()


def test_scvi_get_latent_rep_backwards_compat(n_latent: int = 5):
    adata = synthetic_iid()
    SCVI.setup_anndata(
        adata,
        batch_key="batch",
        labels_key="labels",
    )
    model = SCVI(adata, n_latent=n_latent)
    model.train(1, check_val_every_n_epoch=1, train_size=0.5)
    vae = model.module
    vae_mock = mock.Mock(wraps=model.module)

    def old_inference(*args, **kwargs):
        inf_outs = vae.inference(*args, **kwargs)
        qz = inf_outs.pop("qz")
        inf_outs["qz_m"], inf_outs["qz_v"] = qz.loc, qz.scale**2
        return inf_outs

    vae_mock.inference.side_effect = old_inference
    model.module = vae_mock

    model.get_latent_representation()


def test_scvi_get_feature_corr_backwards_compat(n_latent: int = 5):
    adata = synthetic_iid()
    SCVI.setup_anndata(
        adata,
        batch_key="batch",
        labels_key="labels",
    )
    model = SCVI(adata, n_latent=n_latent)
    model.train(1, check_val_every_n_epoch=1, train_size=0.5)
    vae = model.module
    vae_mock = mock.Mock(wraps=model.module)

    def old_forward(*args, **kwargs):
        inf_outs, gen_outs = vae.forward(*args, **kwargs)
        qz = inf_outs.pop("qz")
        inf_outs["qz_m"], inf_outs["qz_v"] = qz.loc, qz.scale**2
        px = gen_outs.pop("px")
        gen_outs["px_scale"], gen_outs["px_r"] = px.scale, px.theta
        return inf_outs, gen_outs

    vae_mock.forward.side_effect = old_forward
    vae_mock.generative.__signature__ = inspect.signature(
        vae.generative
    )  # Necessary to pass transform_batch check.
    model.module = vae_mock

    model.get_feature_correlation_matrix()


def test_scvi_sparse(n_latent: int = 5):
    adata = synthetic_iid()
    adata.X = csr_matrix(adata.X)
    SCVI.setup_anndata(adata)
    model = SCVI(adata, n_latent=n_latent)
    model.train(1, train_size=0.5)
    assert model.is_trained is True
    z = model.get_latent_representation()
    assert z.shape == (adata.shape[0], n_latent)
    model.get_elbo()
    model.get_marginal_ll(n_mc_samples=3)
    model.get_reconstruction_error()
    model.get_normalized_expression()
    model.differential_expression(groupby="labels", group1="label_1")


def test_setting_adata_attr(n_latent: int = 5):
    adata = synthetic_iid()
    SCVI.setup_anndata(adata, batch_key="batch")
    model = SCVI(adata, n_latent=n_latent)
    model.train(1, train_size=0.5)

    adata2 = synthetic_iid()
    model.adata = adata2

    with pytest.raises(AssertionError):
        rep = model.get_latent_representation(adata)
        rep2 = model.get_latent_representation()
        np.testing.assert_array_equal(rep, rep2)

    orig_manager = model.get_anndata_manager(adata)
    assert model.registry_ is not orig_manager.registry
    assert model.summary_stats is not orig_manager.summary_stats

    adata3 = synthetic_iid()
    del adata3.obs["batch"]
    # validation catches no batch
    with pytest.raises(KeyError):
        model.adata = adata3
        model.get_latent_representation()


def assert_dict_is_subset(d1, d2):
    if not isinstance(d1, dict):
        raise AssertionError(f"{d1} is not a dictionary.")
    elif not isinstance(d2, dict):
        raise AssertionError(f"{d2} is not a dictionary.")

    for k, v in d1.items():
        if k not in d2:
            raise AssertionError(f"{k} missing from {d2}.")
        v2 = d2[k]
        if isinstance(v, dict):
            assert_dict_is_subset(v, v2)
        elif isinstance(v, np.ndarray):
            np.testing.assert_array_equal(v, v2)
        elif v != v2:
            raise AssertionError(f"Mismatch between {v} and {v2}.")


def test_new_setup_compat():
    adata = synthetic_iid()
    adata.obs["cat1"] = np.random.randint(0, 5, size=(adata.shape[0],))
    adata.obs["cat2"] = np.random.randint(0, 5, size=(adata.shape[0],))
    adata.obs["cont1"] = np.random.normal(size=(adata.shape[0],))
    adata.obs["cont2"] = np.random.normal(size=(adata.shape[0],))
    # Handle edge case where registry_key != obs_key.
    adata.obs.rename(
        columns={"batch": "testbatch", "labels": "testlabels"}, inplace=True
    )
    adata2 = adata.copy()

    SCVI.setup_anndata(
        adata,
        batch_key="testbatch",
        labels_key="testlabels",
        categorical_covariate_keys=["cat1", "cat2"],
        continuous_covariate_keys=["cont1", "cont2"],
    )
    model = SCVI(adata)
    adata_manager = model.adata_manager
    model.view_anndata_setup(hide_state_registries=True)

    field_registries = adata_manager.registry[_constants._FIELD_REGISTRIES_KEY]
    field_registries_legacy_subset = {
        k: v for k, v in field_registries.items() if k in LEGACY_REGISTRY_KEYS
    }

    # Backwards compatibility test.
    registry = registry_from_setup_dict(SCVI, LEGACY_SETUP_DICT)
    assert_dict_is_subset(
        registry[_constants._FIELD_REGISTRIES_KEY],
        field_registries_legacy_subset,
    )

    # Test transfer.
    adata2_manager = adata_manager.transfer_fields(adata2)
    np.testing.assert_equal(
        field_registries,
        adata2_manager.registry[_constants._FIELD_REGISTRIES_KEY],
    )


@pytest.mark.internet
def test_backwards_compatible_loading(save_path):
    def download_080_models(save_path):
        file_path = (
            "https://github.com/yoseflab/scVI-data/raw/master/testing_models.tar.gz"
        )
        save_fn = "testing_models.tar.gz"
        _download(file_path, save_path, save_fn)
        saved_file_path = os.path.join(save_path, save_fn)
        tar = tarfile.open(saved_file_path, "r:gz")
        tar.extractall(path=save_path)
        tar.close()

    download_080_models(save_path)
    pretrained_scvi_path = os.path.join(save_path, "testing_models/080_scvi")
    pretrained_scvi_updated_path = os.path.join(
        save_path, "testing_models/080_scvi_updated"
    )
    a = synthetic_iid()
    # Fail legacy load.
    with pytest.raises(ValueError):
        m = SCVI.load(pretrained_scvi_path, adata=a)
    SCVI.convert_legacy_save(pretrained_scvi_path, pretrained_scvi_updated_path)
    m = SCVI.load(pretrained_scvi_updated_path, adata=a)
    m.train(1)


@pytest.mark.internet
def test_backup_url(save_path):
    backup_path = "https://github.com/yoseflab/scVI-data/raw/master/testing_models_0150"
    a = synthetic_iid()
    a.obs["cat1"] = np.random.randint(0, 5, size=(a.shape[0],))
    a.obs["cat2"] = np.random.randint(0, 5, size=(a.shape[0],))
    a.obs["cont1"] = np.random.normal(size=(a.shape[0],))
    a.obs["cont2"] = np.random.normal(size=(a.shape[0],))

    # SCVI
    pretrained_scvi_path = os.path.join(save_path, "testing_models/0150_scvi")
    scvi_backup_url = os.path.join(backup_path, "0150_scvi/model.pt")
    m = SCVI.load(pretrained_scvi_path, adata=a, backup_url=scvi_backup_url)
    m.train(1)


def test_backed_anndata_scvi(save_path):
    adata = synthetic_iid()
    path = os.path.join(save_path, "test_data.h5ad")
    adata.write_h5ad(path)
    adata = anndata.read_h5ad(path, backed="r+")
    SCVI.setup_anndata(adata, batch_key="batch")

    model = SCVI(adata, n_latent=5)
    model.train(1, train_size=0.5)
    assert model.is_trained is True
    z = model.get_latent_representation()
    assert z.shape == (adata.shape[0], 5)
    model.get_elbo()


def test_multiple_covariates_scvi():
    adata = synthetic_iid()
    adata.obs["cont1"] = np.random.normal(size=(adata.shape[0],))
    adata.obs["cont2"] = np.random.normal(size=(adata.shape[0],))
    adata.obs["cat1"] = np.random.randint(0, 5, size=(adata.shape[0],))
    adata.obs["cat2"] = np.random.randint(0, 5, size=(adata.shape[0],))

    SCVI.setup_anndata(
        adata,
        batch_key="batch",
        labels_key="labels",
        continuous_covariate_keys=["cont1", "cont2"],
        categorical_covariate_keys=["cat1", "cat2"],
    )
    m = SCVI(adata)
    m.train(1)
    m.get_latent_representation()
    m.get_elbo()
    m.get_marginal_ll(n_mc_samples=3)
    m.get_reconstruction_error()
    m.get_normalized_expression(n_samples=1)
    m.get_normalized_expression(n_samples=2)


def test_multiple_encoded_covariates_scvi():
    adata = synthetic_iid()
    adata.obs["cont1"] = np.random.normal(size=(adata.shape[0],))
    adata.obs["cont2"] = np.random.normal(size=(adata.shape[0],))
    adata.obs["cat1"] = np.random.randint(0, 5, size=(adata.shape[0],))
    adata.obs["cat2"] = np.random.randint(0, 5, size=(adata.shape[0],))

    SCVI.setup_anndata(
        adata,
        batch_key="batch",
        labels_key="labels",
        continuous_covariate_keys=["cont1", "cont2"],
        categorical_covariate_keys=["cat1", "cat2"],
    )
    m = SCVI(adata, encode_covariates=True)
    m.train(1)


def test_early_stopping():
    n_epochs = 100

    adata = synthetic_iid()
    SCVI.setup_anndata(
        adata,
        batch_key="batch",
        labels_key="labels",
    )
    model = SCVI(adata)
    model.train(n_epochs, early_stopping=True, plan_kwargs={"lr": 0})
    assert len(model.history["elbo_train"]) < n_epochs


def test_de_features():
    adata = synthetic_iid()
    SCVI.setup_anndata(
        adata,
        batch_key="batch",
        labels_key="labels",
    )
    model = SCVI(adata)
    model.train(1)

    model.differential_expression(
        groupby="labels",
        pseudocounts=1e-4,
    )
    model.differential_expression(
        groupby="labels",
        weights="importance",
    )
    model.differential_expression(
        groupby="labels",
        delta=0.5,
        weights="importance",
    )


def test_scarches_data_prep(save_path):
    n_latent = 5
    adata1 = synthetic_iid()
    SCVI.setup_anndata(adata1, batch_key="batch", labels_key="labels")
    model = SCVI(adata1, n_latent=n_latent)
    model.train(1, check_val_every_n_epoch=1)
    dir_path = os.path.join(save_path, "saved_model/")
    model.save(dir_path, overwrite=True)

    # adata2 has more genes and a perfect subset of adata1
    adata2 = synthetic_iid(n_genes=110)
    adata2.obs["batch"] = adata2.obs.batch.cat.rename_categories(["batch_2", "batch_3"])
    SCVI.prepare_query_anndata(adata2, dir_path)
    SCVI.load_query_data(adata2, dir_path)

    adata3 = SCVI.prepare_query_anndata(adata2, dir_path, inplace=False)
    SCVI.load_query_data(adata3, dir_path)

    # adata4 has more genes and missing 10 genes from adata1
    adata4 = synthetic_iid(n_genes=110)
    new_var_names_init = [f"Random {i}" for i in range(10)]
    new_var_names = new_var_names_init + adata4.var_names[10:].to_list()
    adata4.var_names = new_var_names

    SCVI.prepare_query_anndata(adata4, dir_path)
    # should be padded 0s
    assert np.sum(adata4[:, adata4.var_names[:10]].X) == 0
    np.testing.assert_equal(
        adata4.var_names[:10].to_numpy(), adata1.var_names[:10].to_numpy()
    )
    SCVI.load_query_data(adata4, dir_path)

    adata5 = SCVI.prepare_query_anndata(adata4, dir_path, inplace=False)
    SCVI.load_query_data(adata5, dir_path)


def test_scarches_data_prep_layer(save_path):
    n_latent = 5
    adata1 = synthetic_iid()
    adata1.layers["counts"] = adata1.X.copy()
    SCVI.setup_anndata(adata1, layer="counts", batch_key="batch", labels_key="labels")
    model = SCVI(adata1, n_latent=n_latent)
    model.train(1, check_val_every_n_epoch=1)
    dir_path = os.path.join(save_path, "saved_model/")
    model.save(dir_path, overwrite=True)

    # adata4 has more genes and missing 10 genes from adata1
    adata4 = synthetic_iid(n_genes=110)
    adata4.layers["counts"] = adata4.X.copy()
    new_var_names_init = [f"Random {i}" for i in range(10)]
    new_var_names = new_var_names_init + adata4.var_names[10:].to_list()
    adata4.var_names = new_var_names

    SCVI.prepare_query_anndata(adata4, dir_path)
    # should be padded 0s
    assert np.sum(adata4[:, adata4.var_names[:10]].layers["counts"]) == 0
    np.testing.assert_equal(
        adata4.var_names[:10].to_numpy(), adata1.var_names[:10].to_numpy()
    )
    SCVI.load_query_data(adata4, dir_path)


def single_pass_for_online_update(model):
    dl = model._make_data_loader(model.adata, indices=range(0, 10))
    for tensors in dl:
        _, _, scvi_loss = model.module(tensors)
    scvi_loss.loss.backward()


def test_scvi_online_update(save_path):
    n_latent = 5
    adata1 = synthetic_iid()
    SCVI.setup_anndata(adata1, batch_key="batch", labels_key="labels")
    model = SCVI(adata1, n_latent=n_latent)
    model.train(1, check_val_every_n_epoch=1)
    dir_path = os.path.join(save_path, "saved_model/")
    model.save(dir_path, overwrite=True)

    # also test subset var option
    adata2 = synthetic_iid(n_genes=110)
    adata2.obs["batch"] = adata2.obs.batch.cat.rename_categories(["batch_2", "batch_3"])

    model2 = SCVI.load_query_data(adata2, dir_path, inplace_subset_query_vars=True)
    model2.train(max_epochs=1, plan_kwargs={"weight_decay": 0.0})
    model2.get_latent_representation()

    # encoder linear layer equal
    one = (
        model.module.z_encoder.encoder.fc_layers[0][0]
        .weight.detach()
        .cpu()
        .numpy()[:, : adata1.shape[1]]
    )
    two = (
        model2.module.z_encoder.encoder.fc_layers[0][0]
        .weight.detach()
        .cpu()
        .numpy()[:, : adata1.shape[1]]
    )
    np.testing.assert_equal(one, two)
    single_pass_for_online_update(model2)
    assert (
        np.sum(
            model2.module.z_encoder.encoder.fc_layers[0][0]
            .weight.grad.cpu()
            .numpy()[:, : adata1.shape[1]]
        )
        == 0
    )
    # dispersion
    assert model2.module.px_r.requires_grad is False
    # library encoder linear layer
    assert model2.module.l_encoder.encoder.fc_layers[0][0].weight.requires_grad is True
    # 5 for n_latent, 4 for batches
    assert model2.module.decoder.px_decoder.fc_layers[0][0].weight.shape[1] == 9

    # test options
    adata1 = synthetic_iid()
    SCVI.setup_anndata(adata1, batch_key="batch", labels_key="labels")
    model = SCVI(
        adata1,
        n_latent=n_latent,
        n_layers=2,
        encode_covariates=True,
        use_batch_norm="encoder",
        use_layer_norm="none",
    )
    model.train(1, check_val_every_n_epoch=1)
    dir_path = os.path.join(save_path, "saved_model/")
    model.save(dir_path, overwrite=True)

    adata2 = synthetic_iid()
    adata2.obs["batch"] = adata2.obs.batch.cat.rename_categories(["batch_2", "batch_3"])

    model2 = SCVI.load_query_data(adata2, dir_path, freeze_expression=True)
    model2.train(max_epochs=1, plan_kwargs={"weight_decay": 0.0})
    # deactivate no grad decorator
    model2.get_latent_representation()
    # pytorch lightning zeros the grad, so this will get a grad to inspect
    single_pass_for_online_update(model2)
    grad = model2.module.z_encoder.encoder.fc_layers[0][0].weight.grad.cpu().numpy()
    # expression part has zero grad
    assert np.sum(grad[:, :-4]) == 0
    # categorical part has non-zero grad
    assert np.sum(grad[:, -4:]) != 0
    grad = model2.module.decoder.px_decoder.fc_layers[0][0].weight.grad.cpu().numpy()
    # linear layer weight in decoder layer has non-zero grad
    assert np.sum(grad[:, :-4]) == 0

    # do not freeze expression
    model3 = SCVI.load_query_data(
        adata2,
        dir_path,
        freeze_expression=False,
        freeze_batchnorm_encoder=True,
        freeze_decoder_first_layer=False,
    )
    model3.train(max_epochs=1)
    model3.get_latent_representation()
    assert model3.module.z_encoder.encoder.fc_layers[0][1].momentum == 0
    # batch norm weight in encoder layer
    assert model3.module.z_encoder.encoder.fc_layers[0][1].weight.requires_grad is False
    single_pass_for_online_update(model3)
    grad = model3.module.z_encoder.encoder.fc_layers[0][0].weight.grad.cpu().numpy()
    # linear layer weight in encoder layer has non-zero grad
    assert np.sum(grad[:, :-4]) != 0
    grad = model3.module.decoder.px_decoder.fc_layers[0][0].weight.grad.cpu().numpy()
    # linear layer weight in decoder layer has non-zero grad
    assert np.sum(grad[:, :-4]) != 0

    # do not freeze batchnorm
    model3 = SCVI.load_query_data(adata2, dir_path, freeze_batchnorm_encoder=False)
    model3.train(max_epochs=1)
    model3.get_latent_representation()


def test_scvi_library_size_update(save_path):
    n_latent = 5
    adata1 = synthetic_iid()
    SCVI.setup_anndata(adata1, batch_key="batch", labels_key="labels")
    model = SCVI(adata1, n_latent=n_latent, use_observed_lib_size=False)

    assert (
        getattr(model.module, "library_log_means", None) is not None
        and model.module.library_log_means.shape == (1, 2)
        and model.module.library_log_means.count_nonzero().item() == 2
    )
    assert getattr(
        model.module, "library_log_vars", None
    ) is not None and model.module.library_log_vars.shape == (1, 2)

    model.train(1, check_val_every_n_epoch=1)
    dir_path = os.path.join(save_path, "saved_model/")
    model.save(dir_path, overwrite=True)

    # also test subset var option
    adata2 = synthetic_iid(n_genes=110)
    adata2.obs["batch"] = adata2.obs.batch.cat.rename_categories(["batch_2", "batch_3"])

    model2 = SCVI.load_query_data(adata2, dir_path, inplace_subset_query_vars=True)
    assert (
        getattr(model2.module, "library_log_means", None) is not None
        and model2.module.library_log_means.shape == (1, 4)
        and model2.module.library_log_means[:, :2].equal(model.module.library_log_means)
        and model2.module.library_log_means.count_nonzero().item() == 4
    )
    assert (
        getattr(model2.module, "library_log_vars", None) is not None
        and model2.module.library_log_vars.shape == (1, 4)
        and model2.module.library_log_vars[:, :2].equal(model.module.library_log_vars)
    )


def test_set_seed(n_latent: int = 5, seed: int = 1):
    scvi.settings.seed = seed
    adata = synthetic_iid()
    SCVI.setup_anndata(adata, batch_key="batch", labels_key="labels")
    model1 = SCVI(adata, n_latent=n_latent)
    model1.train(1)
    scvi.settings.seed = seed
    model2 = SCVI(adata, n_latent=n_latent)
    model2.train(1)
    assert torch.equal(
        model1.module.z_encoder.encoder.fc_layers[0][0].weight,
        model2.module.z_encoder.encoder.fc_layers[0][0].weight,
    )
