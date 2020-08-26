import pandas as pd

import pytest
from scvi.dataset import synthetic_iid, transfer_anndata_setup, setup_anndata
from scvi.models import SCVI, SCANVI, GIMVI, TOTALVI, LinearSCVI, AUTOZI


def test_SCVI():
    adata = synthetic_iid()
    model = SCVI(adata, n_latent=10)
    model.train(1)
    assert model.is_trained is True
    z = model.get_latent_representation()
    assert z.shape == (adata.shape[0], 10)
    model.get_elbo()
    model.get_marginal_ll()
    model.get_reconstruction_error()
    model.get_normalized_expression()

    adata2 = synthetic_iid()
    model.get_elbo(adata2)
    model.get_marginal_ll(adata2)
    model.get_reconstruction_error(adata2)
    model.get_normalized_expression(adata2)

    # test transfer_anndata_setup
    adata2 = synthetic_iid(run_setup_anndata=False)
    transfer_anndata_setup(adata, adata2)
    model.get_elbo(adata2)

    # test automatic transfer_anndata_setup + on a view
    adata = synthetic_iid()
    model = SCVI(adata)
    adata2 = synthetic_iid(run_setup_anndata=False)
    model.get_elbo(adata2[:10])

    # test that we catch incorrect mappings
    adata = synthetic_iid()
    adata2 = synthetic_iid(run_setup_anndata=False)
    transfer_anndata_setup(adata, adata2)
    adata2.uns["_scvi"]["categorical_mappings"]["_scvi_labels"]["mapping"] = pd.Index(
        data=["undefined_1", "undefined_0", "undefined_2"]
    )
    with pytest.raises(AssertionError):
        model.get_elbo(adata2)


def test_SCANVI():
    adata = synthetic_iid()
    model = SCANVI(adata, "undefined_0", n_latent=10)
    model.train(1)


def test_LinearSCVI():
    # test using raw
    adata = synthetic_iid()
    adata.raw = adata
    adata = adata[:, :10].copy()
    setup_anndata(adata, use_raw=True)
    model = LinearSCVI(adata, n_latent=10)
    model.train(1)
    loadings = model.get_loadings()
    pd.testing.assert_index_equal(loadings.index, adata.raw.var_names)


def test_GIMVI():
    adata = synthetic_iid()
    adata2 = synthetic_iid()
    model = GIMVI(adata, adata2, n_latent=10)
    model.train(1)


def test_AUTOZI():
    data = synthetic_iid(n_batches=1)

    for disp_zi in ["gene", "gene-label"]:
        autozivae = AUTOZI(
            data,
            dispersion=disp_zi,
            zero_inflation=disp_zi,
        )
        autozivae.train(1, lr=1e-2)
        autozivae.get_elbo(indices=autozivae.test_indices)
        autozivae.get_reconstruction_error(indices=autozivae.test_indices)
        autozivae.get_marginal_ll(indices=autozivae.test_indices)
        autozivae.get_alphas_betas()


def test_TOTALVI():
    adata = synthetic_iid()
    n_obs = adata.n_obs
    n_vars = adata.n_vars
    n_proteins = adata.obsm["protein_expression"].shape[1]
    n_latent = 10

    model = TOTALVI(adata, n_latent=n_latent)
    model.train(1)
    assert model.is_trained is True
    z = model.get_latent_representation()
    assert z.shape == (n_obs, n_latent)
    model.get_elbo()
    model.get_marginal_ll()
    model.get_reconstruction_error()
    model.get_normalized_expression()
    model.get_latent_library_size()
    model.get_protein_foreground_probability()
    post_pred = model.posterior_predictive_sample(n_samples=2)
    assert post_pred.shape == (n_obs, n_vars + n_proteins, 2)
    post_pred = model.posterior_predictive_sample(n_samples=1)
    assert post_pred.shape == (n_obs, n_vars + n_proteins)
    feature_correlation_matrix1 = model.get_feature_correlation_matrix(
        correlation_type="spearman"
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
    # model.get_likelihood_parameters()

    model.get_elbo(indices=model.test_indices)
    model.get_marginal_ll(indices=model.test_indices)
    model.get_reconstruction_error(indices=model.test_indices)

    adata2 = synthetic_iid()
    norm_exp = model.get_normalized_expression(adata2, indices=[1, 2, 3])
    assert norm_exp[0].shape == (3, adata2.n_vars)
    assert norm_exp[1].shape == (3, adata2.obsm["protein_expression"].shape[1])

    latent_lib_size = model.get_latent_library_size(adata2, indices=[1, 2, 3])
    assert latent_lib_size.shape == (3, 1)

    pro_foreground_prob = model.get_protein_foreground_probability(
        adata2, indices=[1, 2, 3], protein_list=["1", "2"]
    )
    assert pro_foreground_prob.shape == (3, 2)
    model.posterior_predictive_sample(adata2)
    model.get_feature_correlation_matrix(adata2)
    # model.get_likelihood_parameters(adata2)

    # test transfer_anndata_setup + view
    adata2 = synthetic_iid(run_setup_anndata=False)
    transfer_anndata_setup(adata, adata2)
    model.get_elbo(adata2[:10])

    # test automatic transfer_anndata_setup
    adata = synthetic_iid()
    model = SCVI(adata)
    adata2 = synthetic_iid(run_setup_anndata=False)
    model.get_elbo(adata2)

    # test that we catch incorrect mappings
    adata = synthetic_iid()
    adata2 = synthetic_iid(run_setup_anndata=False)
    transfer_anndata_setup(adata, adata2)
    adata2.uns["_scvi"]["categorical_mappings"]["_scvi_labels"]["mapping"] = pd.Index(
        data=["undefined_1", "undefined_0", "undefined_2"]
    )
    with pytest.raises(AssertionError):
        model.get_elbo(adata2)

    # test that we catch missing proteins
    adata2 = synthetic_iid(run_setup_anndata=False)
    del adata2.obsm["protein_expression"]
    with pytest.raises(KeyError):
        model.get_elbo(adata2)
