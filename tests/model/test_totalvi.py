import os
import pickle

import numpy as np
import pytest
import torch
from anndata import AnnData
from mudata import MuData

from scvi import REGISTRY_KEYS
from scvi.data import pbmcs_10x_cite_seq, synthetic_iid
from scvi.model import TOTALVI
from scvi.utils import attrdict


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
        cls.setup_anndata(
            adata,
            batch_key="batch",
            protein_expression_obsm_key="protein_expression",
            protein_names_uns_key="protein_names",
        )
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

    test_save_load_model(TOTALVI, adata, save_path, prefix=f"{TOTALVI.__name__}_")


@pytest.mark.internet
def test_backup_url(save_path):
    backup_path = "https://github.com/yoseflab/scVI-data/raw/master/testing_models_0150"
    a = synthetic_iid()
    a.obs["cat1"] = np.random.randint(0, 5, size=(a.shape[0],))
    a.obs["cat2"] = np.random.randint(0, 5, size=(a.shape[0],))
    a.obs["cont1"] = np.random.normal(size=(a.shape[0],))
    a.obs["cont2"] = np.random.normal(size=(a.shape[0],))

    # TOTALVI
    pretrained_totalvi_path = os.path.join(save_path, "testing_models/0150_totalvi")
    totalvi_backup_url = os.path.join(backup_path, "0150_totalvi/model.pt")
    m = TOTALVI.load(pretrained_totalvi_path, adata=a, backup_url=totalvi_backup_url)
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
    adata = pbmcs_10x_cite_seq(
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


def test_multiple_covariates_totalvi():
    adata = synthetic_iid()
    adata.obs["cont1"] = np.random.normal(size=(adata.shape[0],))
    adata.obs["cont2"] = np.random.normal(size=(adata.shape[0],))
    adata.obs["cat1"] = np.random.randint(0, 5, size=(adata.shape[0],))
    adata.obs["cat2"] = np.random.randint(0, 5, size=(adata.shape[0],))

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


def test_multiple_encoded_covariates_totalvi():
    adata = synthetic_iid()
    adata.obs["cont1"] = np.random.normal(size=(adata.shape[0],))
    adata.obs["cont2"] = np.random.normal(size=(adata.shape[0],))
    adata.obs["cat1"] = np.random.randint(0, 5, size=(adata.shape[0],))
    adata.obs["cat2"] = np.random.randint(0, 5, size=(adata.shape[0],))

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


def test_totalvi_mudata():
    adata = synthetic_iid()
    protein_adata = synthetic_iid(n_genes=50)
    mdata = MuData({"rna": adata, "protein": protein_adata})
    TOTALVI.setup_mudata(
        mdata,
        batch_key="batch",
        modalities={"rna_layer": "rna", "batch_key": "rna", "protein_layer": "protein"},
    )

    n_obs = mdata.n_obs
    n_genes = adata.n_vars
    n_proteins = protein_adata.X.shape[1]
    n_latent = 10

    model = TOTALVI(mdata, n_latent=n_latent)
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
    assert post_pred.shape == (n_obs, n_genes + n_proteins, 2)
    post_pred = model.posterior_predictive_sample(n_samples=1)
    assert post_pred.shape == (n_obs, n_genes + n_proteins)
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
        n_genes + n_proteins,
        n_genes + n_proteins,
    )
    assert feature_correlation_matrix2.shape == (
        n_genes + n_proteins,
        n_genes + n_proteins,
    )

    model.get_elbo(indices=model.validation_indices)
    model.get_marginal_ll(indices=model.validation_indices, n_mc_samples=3)
    model.get_reconstruction_error(indices=model.validation_indices)

    adata2 = synthetic_iid()
    protein_adata2 = synthetic_iid(n_genes=50)
    mdata2 = MuData({"rna": adata, "protein": protein_adata})
    TOTALVI.setup_mudata(
        mdata2,
        batch_key="batch",
        modalities={"rna_layer": "rna", "batch_key": "rna", "protein_layer": "protein"},
    )
    norm_exp = model.get_normalized_expression(mdata2, indices=[1, 2, 3])
    assert norm_exp[0].shape == (3, adata2.n_vars)
    assert norm_exp[1].shape == (3, protein_adata2.n_vars)
    norm_exp = model.get_normalized_expression(
        mdata2,
        gene_list=adata2.var_names[:5].to_list(),
        protein_list=protein_adata2.var_names[:3].to_list(),
        transform_batch=["batch_0", "batch_1"],
    )

    latent_lib_size = model.get_latent_library_size(mdata2, indices=[1, 2, 3])
    assert latent_lib_size.shape == (3, 1)

    pro_foreground_prob = model.get_protein_foreground_probability(
        mdata2, indices=[1, 2, 3], protein_list=["1", "2"]
    )
    assert pro_foreground_prob.shape == (3, 2)
    model.posterior_predictive_sample(mdata2)
    model.get_feature_correlation_matrix(mdata2)

    # test transfer_anndata_setup + view
    adata2 = synthetic_iid()
    protein_adata2 = synthetic_iid(n_genes=50)
    mdata2 = MuData({"rna": adata2, "protein": protein_adata2})
    model.get_elbo(mdata2[:10])


def test_totalvi_auto_transfer_mudata():
    # test automatic transfer_fields
    adata = synthetic_iid()
    protein_adata = synthetic_iid(n_genes=50)
    mdata = MuData({"rna": adata, "protein": protein_adata})
    TOTALVI.setup_mudata(
        mdata,
        batch_key="batch",
        modalities={"rna_layer": "rna", "batch_key": "rna", "protein_layer": "protein"},
    )
    model = TOTALVI(mdata)
    adata2 = synthetic_iid()
    protein_adata2 = synthetic_iid(n_genes=50)
    mdata2 = MuData({"rna": adata2, "protein": protein_adata2})
    model.get_elbo(mdata2)


def test_totalvi_incorrect_mapping_mudata():
    # test that we catch incorrect mappings
    adata = synthetic_iid()
    protein_adata = synthetic_iid(n_genes=50)
    mdata = MuData({"rna": adata, "protein": protein_adata})
    TOTALVI.setup_mudata(
        mdata,
        batch_key="batch",
        modalities={"rna_layer": "rna", "batch_key": "rna", "protein_layer": "protein"},
    )
    model = TOTALVI(mdata)
    adata2 = synthetic_iid()
    protein_adata2 = synthetic_iid(n_genes=50)
    mdata2 = MuData({"rna": adata2, "protein": protein_adata2})
    adata2.obs.batch = adata2.obs.batch.cat.rename_categories(["batch_0", "batch_10"])
    with pytest.raises(ValueError):
        model.get_elbo(mdata2)


def test_totalvi_reordered_mapping_mudata():
    # test that same mapping different order is okay
    adata = synthetic_iid()
    protein_adata = synthetic_iid(n_genes=50)
    mdata = MuData({"rna": adata, "protein": protein_adata})
    TOTALVI.setup_mudata(
        mdata,
        batch_key="batch",
        modalities={"rna_layer": "rna", "batch_key": "rna", "protein_layer": "protein"},
    )
    model = TOTALVI(mdata)
    adata2 = synthetic_iid()
    protein_adata2 = synthetic_iid(n_genes=50)
    mdata2 = MuData({"rna": adata2, "protein": protein_adata2})
    adata2.obs.batch = adata2.obs.batch.cat.rename_categories(["batch_1", "batch_0"])
    model.get_elbo(mdata2)


def test_totalvi_missing_proteins(save_path):
    # test with missing proteins
    adata = pbmcs_10x_cite_seq(
        save_path=save_path,
        protein_join="outer",
    )
    protein_adata = AnnData(adata.obsm["protein_expression"])
    mdata = MuData({"rna": adata, "protein": protein_adata})
    TOTALVI.setup_mudata(
        mdata,
        batch_key="batch",
        modalities={"rna_layer": "rna", "batch_key": "rna", "protein_layer": "protein"},
    )
    model = TOTALVI(mdata)
    assert model.module.protein_batch_mask is not None
    model.train(1, train_size=0.5)

    model = TOTALVI(mdata, override_missing_proteins=True)
    assert model.module.protein_batch_mask is None
    model.train(1, train_size=0.5)


def test_totalvi_model_library_size_mudata():
    adata = synthetic_iid()
    protein_adata = synthetic_iid(n_genes=50)
    mdata = MuData({"rna": adata, "protein": protein_adata})
    TOTALVI.setup_mudata(
        mdata,
        batch_key="batch",
        modalities={"rna_layer": "rna", "batch_key": "rna", "protein_layer": "protein"},
    )

    n_latent = 10
    model = TOTALVI(mdata, n_latent=n_latent, use_observed_lib_size=False)
    assert hasattr(model.module, "library_log_means") and hasattr(
        model.module, "library_log_vars"
    )
    model.train(1, train_size=0.5)
    assert model.is_trained is True
    model.get_elbo()
    model.get_marginal_ll(n_mc_samples=3)
    model.get_latent_library_size()


def test_totalvi_size_factor_mudata():
    adata = synthetic_iid()
    adata.obs["size_factor"] = np.random.randint(1, 5, size=(adata.shape[0],))
    protein_adata = synthetic_iid(n_genes=50)
    mdata = MuData({"rna": adata, "protein": protein_adata})
    TOTALVI.setup_mudata(
        mdata,
        batch_key="batch",
        size_factor_key="size_factor",
        modalities={
            "rna_layer": "rna",
            "batch_key": "rna",
            "protein_layer": "protein",
            "size_factor_key": "rna",
        },
    )

    n_latent = 10

    # Test size_factor_key overrides use_observed_lib_size.
    model = TOTALVI(mdata, n_latent=n_latent, use_observed_lib_size=False)
    assert not hasattr(model.module, "library_log_means") and not hasattr(
        model.module, "library_log_vars"
    )
    assert model.module.use_size_factor_key
    model.train(1, train_size=0.5)

    model = TOTALVI(mdata, n_latent=n_latent, use_observed_lib_size=True)
    assert not hasattr(model.module, "library_log_means") and not hasattr(
        model.module, "library_log_vars"
    )
    assert model.module.use_size_factor_key
    model.train(1, train_size=0.5)


def test_totalvi_saving_and_loading_mudata(save_path):
    adata = synthetic_iid()
    protein_adata = synthetic_iid(n_genes=50)
    mdata = MuData({"rna": adata, "protein": protein_adata})
    TOTALVI.setup_mudata(
        mdata,
        batch_key="batch",
        modalities={"rna_layer": "rna", "batch_key": "rna", "protein_layer": "protein"},
    )
    model = TOTALVI(mdata)
    model.train(1, train_size=0.2)
    z1 = model.get_latent_representation(mdata)
    test_idx1 = model.validation_indices

    model.save(save_path, overwrite=True, save_anndata=True)
    model.view_setup_args(save_path)

    model = TOTALVI.load(save_path)
    model.get_latent_representation()

    # Load with mismatched genes.
    tmp_adata = synthetic_iid(
        n_genes=200,
    )
    tmp_protein_adata = synthetic_iid(n_genes=50)
    tmp_mdata = MuData({"rna": tmp_adata, "protein": tmp_protein_adata})
    with pytest.raises(ValueError):
        TOTALVI.load(save_path, adata=tmp_mdata)

    # Load with different batches.
    tmp_adata = synthetic_iid()
    tmp_adata.obs["batch"] = tmp_adata.obs["batch"].cat.rename_categories(
        ["batch_2", "batch_3"]
    )
    tmp_protein_adata = synthetic_iid(n_genes=50)
    tmp_mdata = MuData({"rna": tmp_adata, "protein": tmp_protein_adata})
    with pytest.raises(ValueError):
        TOTALVI.load(save_path, adata=tmp_mdata)

    model = TOTALVI.load(save_path, adata=mdata)
    assert REGISTRY_KEYS.BATCH_KEY in model.adata_manager.data_registry
    assert model.adata_manager.data_registry.batch == attrdict(
        {"mod_key": "rna", "attr_name": "obs", "attr_key": "_scvi_batch"}
    )

    z2 = model.get_latent_representation()
    test_idx2 = model.validation_indices
    np.testing.assert_array_equal(z1, z2)
    np.testing.assert_array_equal(test_idx1, test_idx2)
    assert model.is_trained is True

    save_path = os.path.join(save_path, "tmp")

    adata2 = synthetic_iid()
    protein_adata2 = synthetic_iid(n_genes=50)
    mdata2 = MuData({"rna": adata2, "protein": protein_adata2})
    TOTALVI.setup_mudata(
        mdata2,
        batch_key="batch",
        modalities={"rna_layer": "rna", "batch_key": "rna", "protein_layer": "protein"},
    )


def test_totalvi_online_update(save_path):
    # basic case
    n_latent = 5
    adata1 = synthetic_iid()
    TOTALVI.setup_anndata(
        adata1,
        batch_key="batch",
        protein_expression_obsm_key="protein_expression",
        protein_names_uns_key="protein_names",
    )
    model = TOTALVI(adata1, n_latent=n_latent, use_batch_norm="decoder")
    model.train(1, check_val_every_n_epoch=1)
    dir_path = os.path.join(save_path, "saved_model/")
    model.save(dir_path, overwrite=True)

    adata2 = synthetic_iid()
    adata2.obs["batch"] = adata2.obs.batch.cat.rename_categories(["batch_2", "batch_3"])

    model2 = TOTALVI.load_query_data(adata2, dir_path)
    assert model2.module.background_pro_alpha.requires_grad is True
    model2.train(max_epochs=1)
    model2.get_latent_representation()

    # batch 3 has no proteins
    adata2 = synthetic_iid()
    adata2.obs["batch"] = adata2.obs.batch.cat.rename_categories(["batch_2", "batch_3"])
    adata2.obsm["protein_expression"][adata2.obs.batch == "batch_3"] = 0

    # load from model in memory
    model3 = TOTALVI.load_query_data(adata2, model)
    model3.train(max_epochs=1)
    model3.get_latent_representation()
