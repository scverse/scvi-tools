import os
import pickle
import tarfile

import numpy as np
import pytest
import torch

import scvi
from scvi.data import synthetic_iid
from scvi.data._compat import LEGACY_REGISTRY_KEY_MAP
from scvi.data._download import _download
from scvi.model import (
    AUTOZI,
    PEAKVI,
    SCANVI,
    SCVI,
    TOTALVI,
    LinearSCVI,
)
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
        if cls is TOTALVI:
            cls.setup_anndata(
                adata,
                batch_key="batch",
                protein_expression_obsm_key="protein_expression",
                protein_names_uns_key="protein_names",
            )
        else:
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

    for cls in [SCVI, LinearSCVI, TOTALVI, PEAKVI]:
        test_save_load_model(cls, adata, save_path, prefix=f"{cls.__name__}_")

    # AUTOZI
    prefix = "AUTOZI_"
    AUTOZI.setup_anndata(adata, batch_key="batch", labels_key="labels")
    model = AUTOZI(adata, latent_distribution="normal")
    model.train(1, train_size=0.5)
    ab1 = model.get_alphas_betas()
    model.save(save_path, overwrite=True, save_anndata=True, prefix=prefix)
    model.view_setup_args(save_path, prefix=prefix)
    model = AUTOZI.load(save_path, prefix=prefix)
    model.get_latent_representation()
    tmp_adata = scvi.data.synthetic_iid(n_genes=200)
    with pytest.raises(ValueError):
        AUTOZI.load(save_path, adata=tmp_adata, prefix=prefix)
    model = AUTOZI.load(save_path, adata=adata, prefix=prefix)
    assert "batch" in model.adata_manager.data_registry
    assert model.adata_manager.data_registry.batch == attrdict(
        {"attr_name": "obs", "attr_key": "_scvi_batch"}
    )

    ab2 = model.get_alphas_betas()
    np.testing.assert_array_equal(ab1["alpha_posterior"], ab2["alpha_posterior"])
    np.testing.assert_array_equal(ab1["beta_posterior"], ab2["beta_posterior"])
    assert model.is_trained is True

    # Test legacy loading
    legacy_save_path = os.path.join(save_path, "legacy/")
    legacy_save(
        model, legacy_save_path, overwrite=True, save_anndata=True, prefix=prefix
    )
    with pytest.raises(ValueError):
        AUTOZI.load(legacy_save_path, adata=adata, prefix=prefix)
    AUTOZI.convert_legacy_save(
        legacy_save_path, legacy_save_path, overwrite=True, prefix=prefix
    )
    m = AUTOZI.load(legacy_save_path, adata=adata, prefix=prefix)
    m.train(1)

    # SCANVI
    prefix = "SCANVI_"
    SCANVI.setup_anndata(adata, "labels", "label_0", batch_key="batch")
    model = SCANVI(adata)
    model.train(max_epochs=1, train_size=0.5)
    p1 = model.predict()
    model.save(save_path, overwrite=True, save_anndata=True, prefix=prefix)
    model.view_setup_args(save_path, prefix=prefix)
    model = SCANVI.load(save_path, prefix=prefix)
    model.get_latent_representation()
    tmp_adata = scvi.data.synthetic_iid(n_genes=200)
    with pytest.raises(ValueError):
        SCANVI.load(save_path, adata=tmp_adata, prefix=prefix)
    model = SCANVI.load(save_path, adata=adata, prefix=prefix)
    assert "batch" in model.adata_manager.data_registry
    assert model.adata_manager.data_registry.batch == attrdict(
        {"attr_name": "obs", "attr_key": "_scvi_batch"}
    )

    p2 = model.predict()
    np.testing.assert_array_equal(p1, p2)
    assert model.is_trained is True

    # Test legacy loading
    legacy_save_path = os.path.join(save_path, "legacy/")
    legacy_save(
        model, legacy_save_path, overwrite=True, save_anndata=True, prefix=prefix
    )
    with pytest.raises(ValueError):
        SCANVI.load(legacy_save_path, adata=adata, prefix=prefix)
    SCANVI.convert_legacy_save(
        legacy_save_path, legacy_save_path, overwrite=True, prefix=prefix
    )
    m = SCANVI.load(legacy_save_path, adata=adata, prefix=prefix)
    m.train(1)


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
    a = scvi.data.synthetic_iid()
    # Fail legacy load.
    with pytest.raises(ValueError):
        m = scvi.model.SCVI.load(pretrained_scvi_path, adata=a)
    scvi.model.SCVI.convert_legacy_save(
        pretrained_scvi_path, pretrained_scvi_updated_path
    )
    m = scvi.model.SCVI.load(pretrained_scvi_updated_path, adata=a)
    m.train(1)
    pretrained_totalvi_path = os.path.join(save_path, "testing_models/080_totalvi")
    pretrained_totalvi_updated_path = os.path.join(
        save_path, "testing_models/080_totalvi_updated"
    )
    # Fail legacy load.
    with pytest.raises(ValueError):
        m = scvi.model.TOTALVI.load(pretrained_totalvi_path, adata=a)
    scvi.model.TOTALVI.convert_legacy_save(
        pretrained_totalvi_path, pretrained_totalvi_updated_path
    )
    m = scvi.model.TOTALVI.load(pretrained_totalvi_updated_path, adata=a)
    m.train(1)


@pytest.mark.internet
def test_backup_url(save_path):
    backup_path = "https://github.com/yoseflab/scVI-data/raw/master/testing_models_0150"
    a = scvi.data.synthetic_iid()
    a.obs["cat1"] = np.random.randint(0, 5, size=(a.shape[0],))
    a.obs["cat2"] = np.random.randint(0, 5, size=(a.shape[0],))
    a.obs["cont1"] = np.random.normal(size=(a.shape[0],))
    a.obs["cont2"] = np.random.normal(size=(a.shape[0],))

    # SCVI
    pretrained_scvi_path = os.path.join(save_path, "testing_models/0150_scvi")
    scvi_backup_url = os.path.join(backup_path, "0150_scvi/model.pt")
    m = scvi.model.SCVI.load(pretrained_scvi_path, adata=a, backup_url=scvi_backup_url)
    m.train(1)

    # TOTALVI
    pretrained_totalvi_path = os.path.join(save_path, "testing_models/0150_totalvi")
    totalvi_backup_url = os.path.join(backup_path, "0150_totalvi/model.pt")
    m = scvi.model.TOTALVI.load(
        pretrained_totalvi_path, adata=a, backup_url=totalvi_backup_url
    )
    m.train(1)


def test_totalvi(save_path):
    adata = synthetic_iid()
    TOTALVI.setup_anndata(
        adata,
        batch_key="batch",
        protein_expression_obsm_key="protein_expression",
        protein_names_uns_key="protein_names",
    )

    n_obs = adata.n_obs
    n_vars = adata.n_vars
    n_proteins = adata.obsm["protein_expression"].shape[1]
    n_latent = 10

    model = TOTALVI(adata, n_latent=n_latent)
    model.train(1, train_size=0.5)
    assert model.is_trained is True
    z = model.get_latent_representation()
    assert z.shape == (n_obs, n_latent)
    model.get_elbo()
    model.get_marginal_ll(n_mc_samples=3)
    model.get_reconstruction_error()
    model.get_normalized_expression()
    model.get_normalized_expression(transform_batch=["batch_0", "batch_1"])
    model.get_latent_library_size()
    model.get_protein_foreground_probability()
    model.get_protein_foreground_probability(transform_batch=["batch_0", "batch_1"])
    post_pred = model.posterior_predictive_sample(n_samples=2)
    assert post_pred.shape == (n_obs, n_vars + n_proteins, 2)
    post_pred = model.posterior_predictive_sample(n_samples=1)
    assert post_pred.shape == (n_obs, n_vars + n_proteins)
    feature_correlation_matrix1 = model.get_feature_correlation_matrix(
        correlation_type="spearman"
    )
    feature_correlation_matrix1 = model.get_feature_correlation_matrix(
        correlation_type="spearman", transform_batch=["batch_0", "batch_1"]
    )
    feature_correlation_matrix2 = model.get_feature_correlation_matrix(
        correlation_type="pearson"
    )
    assert feature_correlation_matrix1.shape == (
        n_vars + n_proteins,
        n_vars + n_proteins,
    )
    assert feature_correlation_matrix2.shape == (
        n_vars + n_proteins,
        n_vars + n_proteins,
    )

    model.get_elbo(indices=model.validation_indices)
    model.get_marginal_ll(indices=model.validation_indices, n_mc_samples=3)
    model.get_reconstruction_error(indices=model.validation_indices)

    adata2 = synthetic_iid()
    TOTALVI.setup_anndata(
        adata2,
        batch_key="batch",
        protein_expression_obsm_key="protein_expression",
        protein_names_uns_key="protein_names",
    )
    norm_exp = model.get_normalized_expression(adata2, indices=[1, 2, 3])
    assert norm_exp[0].shape == (3, adata2.n_vars)
    assert norm_exp[1].shape == (3, adata2.obsm["protein_expression"].shape[1])
    norm_exp = model.get_normalized_expression(
        adata2,
        gene_list=adata2.var_names[:5].to_list(),
        protein_list=adata2.uns["protein_names"][:3],
        transform_batch=["batch_0", "batch_1"],
    )

    latent_lib_size = model.get_latent_library_size(adata2, indices=[1, 2, 3])
    assert latent_lib_size.shape == (3, 1)

    pro_foreground_prob = model.get_protein_foreground_probability(
        adata2, indices=[1, 2, 3], protein_list=["1", "2"]
    )
    assert pro_foreground_prob.shape == (3, 2)
    model.posterior_predictive_sample(adata2)
    model.get_feature_correlation_matrix(adata2)

    # test transfer_anndata_setup + view
    adata2 = synthetic_iid()
    model.get_elbo(adata2[:10])

    # test automatic transfer_anndata_setup
    adata = synthetic_iid()
    # no protein names so we test our auto generation
    TOTALVI.setup_anndata(
        adata,
        batch_key="batch",
        protein_expression_obsm_key="protein_expression",
    )
    model = TOTALVI(adata)
    model.train(1, train_size=0.5)
    adata2 = synthetic_iid()
    model.get_elbo(adata2)

    # test that we catch incorrect mappings
    adata2 = synthetic_iid()
    adata2.obs.batch = adata2.obs.batch.cat.rename_categories(["batch_0", "batch_10"])
    with pytest.raises(ValueError):
        model.get_elbo(adata2)

    # test that same mapping different order is okay
    adata2 = synthetic_iid()
    adata2.obs.batch = adata2.obs.batch.cat.rename_categories(["batch_1", "batch_0"])
    model.get_elbo(adata2)  # should automatically transfer setup

    # test that we catch missing proteins
    adata2 = synthetic_iid()
    del adata2.obsm["protein_expression"]
    with pytest.raises(KeyError):
        model.get_elbo(adata2)
    model.differential_expression(groupby="labels", group1="label_1")
    model.differential_expression(groupby="labels", group1="label_1", group2="label_2")
    model.differential_expression(idx1=[0, 1, 2], idx2=[3, 4, 5])
    model.differential_expression(idx1=[0, 1, 2])
    model.differential_expression(groupby="labels")

    # test with missing proteins
    adata = scvi.data.pbmcs_10x_cite_seq(
        save_path=save_path,
        protein_join="outer",
    )
    TOTALVI.setup_anndata(
        adata, batch_key="batch", protein_expression_obsm_key="protein_expression"
    )
    model = TOTALVI(adata)
    assert model.module.protein_batch_mask is not None
    model.train(1, train_size=0.5)

    model = TOTALVI(adata, override_missing_proteins=True)
    assert model.module.protein_batch_mask is None
    model.train(1, train_size=0.5)


def test_totalvi_model_library_size(save_path):
    adata = synthetic_iid()
    TOTALVI.setup_anndata(
        adata,
        batch_key="batch",
        protein_expression_obsm_key="protein_expression",
        protein_names_uns_key="protein_names",
    )
    n_latent = 10

    model = TOTALVI(adata, n_latent=n_latent, use_observed_lib_size=False)
    assert hasattr(model.module, "library_log_means") and hasattr(
        model.module, "library_log_vars"
    )
    model.train(1, train_size=0.5)
    assert model.is_trained is True
    model.get_elbo()
    model.get_marginal_ll(n_mc_samples=3)
    model.get_latent_library_size()


def test_totalvi_size_factor():
    adata = synthetic_iid()
    adata.obs["size_factor"] = np.random.randint(1, 5, size=(adata.shape[0],))
    TOTALVI.setup_anndata(
        adata,
        batch_key="batch",
        protein_expression_obsm_key="protein_expression",
        protein_names_uns_key="protein_names",
        size_factor_key="size_factor",
    )
    n_latent = 10

    # Test size_factor_key overrides use_observed_lib_size.
    model = TOTALVI(adata, n_latent=n_latent, use_observed_lib_size=False)
    assert not hasattr(model.module, "library_log_means") and not hasattr(
        model.module, "library_log_vars"
    )
    assert model.module.use_size_factor_key
    model.train(1, train_size=0.5)

    model = TOTALVI(adata, n_latent=n_latent, use_observed_lib_size=True)
    assert not hasattr(model.module, "library_log_means") and not hasattr(
        model.module, "library_log_vars"
    )
    assert model.module.use_size_factor_key
    model.train(1, train_size=0.5)


def test_multiple_covariates_scvi(save_path):
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

    SCANVI.setup_anndata(
        adata,
        "labels",
        "Unknown",
        batch_key="batch",
        continuous_covariate_keys=["cont1", "cont2"],
        categorical_covariate_keys=["cat1", "cat2"],
    )
    m = SCANVI(adata)
    m.train(1)
    m.get_latent_representation()
    m.get_elbo()
    m.get_marginal_ll(n_mc_samples=3)
    m.get_reconstruction_error()
    m.get_normalized_expression(n_samples=1)
    m.get_normalized_expression(n_samples=2)

    TOTALVI.setup_anndata(
        adata,
        batch_key="batch",
        protein_expression_obsm_key="protein_expression",
        protein_names_uns_key="protein_names",
        continuous_covariate_keys=["cont1", "cont2"],
        categorical_covariate_keys=["cat1", "cat2"],
    )
    m = TOTALVI(adata)
    m.train(1)
    m.get_latent_representation()
    m.get_elbo()
    m.get_marginal_ll(n_mc_samples=3)
    m.get_reconstruction_error()
    m.get_normalized_expression(n_samples=1)
    m.get_normalized_expression(n_samples=2)


def test_multiple_encoded_covariates_scvi(save_path):
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

    SCANVI.setup_anndata(
        adata,
        "labels",
        "Unknown",
        batch_key="batch",
        continuous_covariate_keys=["cont1", "cont2"],
        categorical_covariate_keys=["cat1", "cat2"],
    )
    m = SCANVI(adata, encode_covariates=True)
    m.train(1)

    TOTALVI.setup_anndata(
        adata,
        batch_key="batch",
        protein_expression_obsm_key="protein_expression",
        protein_names_uns_key="protein_names",
        continuous_covariate_keys=["cont1", "cont2"],
        categorical_covariate_keys=["cat1", "cat2"],
    )
    m = TOTALVI(adata, encode_covariates=True)
    m.train(1)
