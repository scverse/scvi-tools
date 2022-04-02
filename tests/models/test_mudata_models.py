import numpy as np
import pytest
from anndata import AnnData
from mudata import MuData

import scvi
from scvi.data import synthetic_iid
from scvi.model import TOTALVI


def test_totalvi(save_path):
    adata = synthetic_iid()
    protein_adata = synthetic_iid(n_genes=50)
    mdata = MuData({"rna": adata, "protein": protein_adata})
    TOTALVI.setup_mudata(
        mdata,
        layer_mod="rna",
        batch_mod="rna",
        batch_key="batch",
        protein_expression_mod="protein",
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
        layer_mod="rna",
        batch_mod="rna",
        batch_key="batch",
        protein_expression_mod="protein",
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


def test_totalvi_auto_transfer(save_path):
    # test automatic transfer_fields
    adata = synthetic_iid()
    protein_adata = synthetic_iid(n_genes=50)
    mdata = MuData({"rna": adata, "protein": protein_adata})
    TOTALVI.setup_mudata(
        mdata,
        layer_mod="rna",
        batch_mod="rna",
        batch_key="batch",
        protein_expression_mod="protein",
    )
    model = TOTALVI(mdata)
    adata2 = synthetic_iid()
    protein_adata2 = synthetic_iid(n_genes=50)
    mdata2 = MuData({"rna": adata2, "protein": protein_adata2})
    model.get_elbo(mdata2)


def test_totalvi_incorrect_mapping(save_path):
    # test that we catch incorrect mappings
    adata = synthetic_iid()
    protein_adata = synthetic_iid(n_genes=50)
    mdata = MuData({"rna": adata, "protein": protein_adata})
    TOTALVI.setup_mudata(
        mdata,
        layer_mod="rna",
        batch_mod="rna",
        batch_key="batch",
        protein_expression_mod="protein",
    )
    model = TOTALVI(mdata)
    adata2 = synthetic_iid()
    protein_adata2 = synthetic_iid(n_genes=50)
    mdata2 = MuData({"rna": adata2, "protein": protein_adata2})
    adata2.obs.batch.cat.rename_categories(["batch_0", "batch_10"], inplace=True)
    with pytest.raises(ValueError):
        model.get_elbo(mdata2)


def test_totalvi_reordered_mapping(save_path):
    # test that same mapping different order is okay
    adata = synthetic_iid()
    protein_adata = synthetic_iid(n_genes=50)
    mdata = MuData({"rna": adata, "protein": protein_adata})
    TOTALVI.setup_mudata(
        mdata,
        layer_mod="rna",
        batch_mod="rna",
        batch_key="batch",
        protein_expression_mod="protein",
    )
    model = TOTALVI(mdata)
    adata2 = synthetic_iid()
    protein_adata2 = synthetic_iid(n_genes=50)
    mdata2 = MuData({"rna": adata2, "protein": protein_adata2})
    adata2.obs.batch.cat.rename_categories(["batch_1", "batch_0"], inplace=True)
    model.get_elbo(mdata2)


def test_totalvi_missing_proteins(save_path):
    # test with missing proteins
    adata = scvi.data.pbmcs_10x_cite_seq(
        save_path=save_path,
        protein_join="outer",
    )
    protein_adata = AnnData(adata.obsm["protein_expression"])
    mdata = MuData({"rna": adata, "protein": protein_adata})
    TOTALVI.setup_mudata(
        mdata,
        layer_mod="rna",
        batch_mod="rna",
        batch_key="batch",
        protein_expression_mod="protein",
    )
    model = TOTALVI(mdata)
    assert model.module.protein_batch_mask is not None
    model.train(1, train_size=0.5)

    model = TOTALVI(mdata, override_missing_proteins=True)
    assert model.module.protein_batch_mask is None
    model.train(1, train_size=0.5)


def test_totalvi_model_library_size(save_path):
    adata = synthetic_iid()
    protein_adata = synthetic_iid(n_genes=50)
    mdata = MuData({"rna": adata, "protein": protein_adata})
    TOTALVI.setup_mudata(
        mdata,
        layer_mod="rna",
        batch_mod="rna",
        batch_key="batch",
        protein_expression_mod="protein",
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


def test_totalvi_size_factor():
    adata = synthetic_iid()
    adata.obs["size_factor"] = np.random.randint(1, 5, size=(adata.shape[0],))
    protein_adata = synthetic_iid(n_genes=50)
    mdata = MuData({"rna": adata, "protein": protein_adata})
    TOTALVI.setup_mudata(
        mdata,
        layer_mod="rna",
        batch_mod="rna",
        batch_key="batch",
        protein_expression_mod="protein",
        size_factor_mod="rna",
        size_factor_key="size_factor",
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
