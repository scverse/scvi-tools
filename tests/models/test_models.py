import pandas as pd
import numpy as np
import pytest

import scvi
from scvi.dataset import synthetic_iid, transfer_anndata_setup, setup_anndata
from scvi.models import SCVI, SCANVI, GIMVI, TOTALVI, LinearSCVI, AUTOZI


def test_SCVI():
    n_latent = 5
    adata = synthetic_iid()
    model = SCVI(adata, n_latent=n_latent)
    model.train(1)
    assert model.is_trained is True
    z = model.get_latent_representation()
    assert z.shape == (adata.shape[0], n_latent)
    model.get_elbo()
    model.get_marginal_ll()
    model.get_reconstruction_error()
    model.get_normalized_expression()

    adata2 = synthetic_iid()
    model.get_elbo(adata2)
    model.get_marginal_ll(adata2)
    model.get_reconstruction_error(adata2)
    latent = model.get_latent_representation(adata2, indices=[1, 2, 3])
    assert latent.shape == (3, n_latent)
    denoised = model.get_normalized_expression(adata2)
    assert denoised.shape == adata.shape

    denoised = model.get_normalized_expression(
        adata2, indices=[1, 2, 3], transform_batch=1
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

    # test mismatched categories raises ValueError
    adata2 = synthetic_iid(run_setup_anndata=False)
    adata2.obs.labels.cat.rename_categories(["a", "b", "c"], inplace=True)
    with pytest.raises(ValueError):
        model.get_elbo(adata2)

    # test differential expression
    model.differential_expression(groupby="labels", group1="undefined_1")
    model.differential_expression(
        groupby="labels", group1="undefined_1", group2="undefined_2", mode="change"
    )


def test_saving_and_loading():
    adata = synthetic_iid()

    for cls in [SCVI, LinearSCVI, TOTALVI]:
        model = cls(adata, latent_distribution="normal")
        model.train(1)
        z1 = model.get_latent_representation(adata)
        test_idx1 = model.test_indices
        model.save("tmp", overwrite=True)
        model = cls.load(adata, "tmp")
        z2 = model.get_latent_representation()
        test_idx2 = model.test_indices
        np.testing.assert_array_equal(z1, z2)
        np.testing.assert_array_equal(test_idx1, test_idx2)
        assert model.is_trained is True

    # AUTOZI
    model = AUTOZI(adata, latent_distribution="normal")
    model.train(1)
    ab1 = model.get_alphas_betas()
    model.save("tmp", overwrite=True)
    model = AUTOZI.load(adata, "tmp")
    ab2 = model.get_alphas_betas()
    np.testing.assert_array_equal(ab1["alpha_posterior"], ab2["alpha_posterior"])
    np.testing.assert_array_equal(ab1["beta_posterior"], ab2["beta_posterior"])
    assert model.is_trained is True

    # SCANVI
    model = SCANVI(adata, "undefined_0")
    model.train(n_epochs_unsupervised=1, n_epochs_semisupervised=1)
    p1 = model.predict()
    model.save("tmp", overwrite=True)
    model = SCANVI.load(adata, "tmp")
    p2 = model.predict()
    np.testing.assert_array_equal(p1, p2)
    assert model.is_trained is True

    # GIMVI
    model = GIMVI(adata, adata)
    model.train(1)
    z1 = model.get_latent_representation([adata])
    z2 = model.get_latent_representation([adata])
    np.testing.assert_array_equal(z1, z2)
    model.save("tmp", overwrite=True)
    model = GIMVI.load(adata, adata, "tmp")
    z2 = model.get_latent_representation([adata])
    np.testing.assert_array_equal(z1, z2)
    assert model.is_trained is True


def test_SCANVI():
    adata = synthetic_iid()
    model = SCANVI(adata, "undefined_0", n_latent=10)
    model.train(1)
    adata2 = synthetic_iid()
    predictions = model.predict(adata2, indices=[1, 2, 3])
    assert len(predictions) == 3
    model.predict(adata2, soft=True)
    model.get_normalized_expression(adata2)
    model.differential_expression(groupby="labels", group1="undefined_1")
    model.differential_expression(
        groupby="labels", group1="undefined_1", group2="undefined_2"
    )


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
    model.differential_expression(groupby="labels", group1="undefined_1")
    model.differential_expression(
        groupby="labels", group1="undefined_1", group2="undefined_2"
    )


def test_GIMVI():
    adata_seq = synthetic_iid()
    adata_spatial = synthetic_iid()
    model = GIMVI(adata_seq, adata_spatial, n_latent=10)
    model.get_latent_representation()
    model.get_imputed_values()
    model.train(1)

    trainer = model.trainer
    results = pd.DataFrame(
        trainer.get_loss_magnitude(),
        index=["reconstruction", "kl_divergence", "discriminator"],
        columns=["Sequencing", "Spatial"],
    )
    results.columns.name = "Dataset"
    results.index.name = "Loss"

    trainer.get_discriminator_confusion()


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


def test_TOTALVI(save_path):
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
    model = TOTALVI(adata)
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
    model.differential_expression(groupby="labels", group1="undefined_1")
    model.differential_expression(
        groupby="labels", group1="undefined_1", group2="undefined_2"
    )

    # test with missing proteins
    adata = scvi.dataset.pbmcs_10x_cite_seq(save_path=save_path, protein_join="outer")
    model = TOTALVI(adata)
    assert model.model.protein_batch_mask is not None
    model.train(2, train_size=0.5)
