import os
import pandas as pd
import numpy as np
import pytest

import scvi
from scvi.data import synthetic_iid, transfer_anndata_setup, setup_anndata
from scvi.model import SCVI, SCANVI, GIMVI, TOTALVI, LinearSCVI, AUTOZI
from scipy.sparse import csr_matrix


def test_scvi():
    n_latent = 5
    adata = synthetic_iid()
    model = SCVI(adata, n_latent=n_latent)
    model.train(1, frequency=1)
    assert model.is_trained is True
    z = model.get_latent_representation()
    assert z.shape == (adata.shape[0], n_latent)
    # len of history should be 2 since metrics is also run once at the very end after training
    assert len(model.history["elbo_train_set"]) == 2
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
    adata2.uns["_scvi"]["categorical_mappings"]["_scvi_labels"]["mapping"] = np.array(
        ["undefined_1", "undefined_0", "undefined_2"]
    )
    with pytest.raises(ValueError):
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
    model.differential_expression(groupby="labels")
    model.differential_expression(idx1=[0, 1, 2], idx2=[3, 4, 5])
    model.differential_expression(idx1=[0, 1, 2])


def test_scvi_sparse():
    n_latent = 5
    adata = synthetic_iid(run_setup_anndata=False)
    adata.X = csr_matrix(adata.X)
    setup_anndata(adata)
    model = SCVI(adata, n_latent=n_latent)
    model.train(1)
    assert model.is_trained is True
    z = model.get_latent_representation()
    assert z.shape == (adata.shape[0], n_latent)
    model.get_elbo()
    model.get_marginal_ll()
    model.get_reconstruction_error()
    model.get_normalized_expression()
    model.differential_expression(groupby="labels", group1="undefined_1")


def test_saving_and_loading(save_path):
    def test_save_load_model(cls, adata, save_path):
        model = cls(adata, latent_distribution="normal")
        model.train(1)
        z1 = model.get_latent_representation(adata)
        test_idx1 = model.test_indices
        model.save(save_path, overwrite=True)
        model = cls.load(adata, save_path)
        z2 = model.get_latent_representation()
        test_idx2 = model.test_indices
        np.testing.assert_array_equal(z1, z2)
        np.testing.assert_array_equal(test_idx1, test_idx2)
        assert model.is_trained is True

    save_path = os.path.join(save_path, "tmp")
    adata = synthetic_iid()

    for cls in [SCVI, LinearSCVI, TOTALVI]:
        print(cls)
        test_save_load_model(cls, adata, save_path)

    # AUTOZI
    model = AUTOZI(adata, latent_distribution="normal")
    model.train(1)
    ab1 = model.get_alphas_betas()
    model.save(save_path, overwrite=True)
    model = AUTOZI.load(adata, save_path)
    ab2 = model.get_alphas_betas()
    np.testing.assert_array_equal(ab1["alpha_posterior"], ab2["alpha_posterior"])
    np.testing.assert_array_equal(ab1["beta_posterior"], ab2["beta_posterior"])
    assert model.is_trained is True

    # SCANVI
    model = SCANVI(adata, "undefined_0")
    model.train(n_epochs_unsupervised=1, n_epochs_semisupervised=1)
    p1 = model.predict()
    model.save(save_path, overwrite=True)
    model = SCANVI.load(adata, save_path)
    p2 = model.predict()
    np.testing.assert_array_equal(p1, p2)
    assert model.is_trained is True

    # GIMVI
    model = GIMVI(adata, adata)
    model.train(1)
    z1 = model.get_latent_representation([adata])
    z2 = model.get_latent_representation([adata])
    np.testing.assert_array_equal(z1, z2)
    model.save(save_path, overwrite=True)
    model = GIMVI.load(adata, adata, save_path)
    z2 = model.get_latent_representation([adata])
    np.testing.assert_array_equal(z1, z2)
    assert model.is_trained is True


def test_scanvi():
    adata = synthetic_iid()
    model = SCANVI(adata, "undefined_0", n_latent=10)
    model.train(1, frequency=1)
    assert len(model.history["unsupervised_trainer_history"]) == 2
    assert len(model.history["semisupervised_trainer_history"]) == 3
    adata2 = synthetic_iid()
    predictions = model.predict(adata2, indices=[1, 2, 3])
    assert len(predictions) == 3
    model.predict()
    model.predict(adata2, soft=True)
    model.predict(adata2, soft=True, indices=[1, 2, 3])
    model.get_normalized_expression(adata2)
    model.differential_expression(groupby="labels", group1="undefined_1")
    model.differential_expression(
        groupby="labels", group1="undefined_1", group2="undefined_2"
    )


def test_linear_scvi():
    # test using raw
    adata = synthetic_iid()
    adata.raw = adata
    adata = adata[:, :10].copy()
    setup_anndata(adata, use_raw=True)
    model = LinearSCVI(adata, n_latent=10)
    model.train(1, frequency=1)
    assert len(model.history["elbo_train_set"]) == 2
    assert len(model.history["elbo_test_set"]) == 2
    loadings = model.get_loadings()
    pd.testing.assert_index_equal(loadings.index, adata.raw.var_names)
    model.differential_expression(groupby="labels", group1="undefined_1")
    model.differential_expression(
        groupby="labels", group1="undefined_1", group2="undefined_2"
    )


def test_gimvi():
    adata_seq = synthetic_iid()
    adata_spatial = synthetic_iid()
    model = GIMVI(adata_seq, adata_spatial, n_latent=10)
    model.get_latent_representation()
    model.get_imputed_values()
    model.train(1, frequency=1, early_stopping_kwargs=None)

    assert len(model.history["elbo_train_0"]) == 2
    assert len(model.history["elbo_train_1"]) == 2
    assert len(model.history["elbo_test_0"]) == 2
    assert len(model.history["elbo_test_1"]) == 2

    trainer = model.trainer
    results = pd.DataFrame(
        trainer.get_loss_magnitude(),
        index=["reconstruction", "kl_divergence", "discriminator"],
        columns=["Sequencing", "Spatial"],
    )
    results.columns.name = "Dataset"
    results.index.name = "Loss"
    trainer.get_discriminator_confusion()
    adata_spatial.var_names += "asdf"
    with pytest.raises(ValueError):
        model = GIMVI(adata_seq, adata_spatial)


def test_autozi():
    data = synthetic_iid(n_batches=1)

    for disp_zi in ["gene", "gene-label"]:
        autozivae = AUTOZI(
            data,
            dispersion=disp_zi,
            zero_inflation=disp_zi,
        )
        autozivae.train(1, lr=1e-2, frequency=1)
        assert len(autozivae.history["elbo_train_set"]) == 2
        assert len(autozivae.history["elbo_test_set"]) == 2
        autozivae.get_elbo(indices=autozivae.test_indices)
        autozivae.get_reconstruction_error(indices=autozivae.test_indices)
        autozivae.get_marginal_ll(indices=autozivae.test_indices)
        autozivae.get_alphas_betas()


def test_totalvi(save_path):
    adata = synthetic_iid()
    n_obs = adata.n_obs
    n_vars = adata.n_vars
    n_proteins = adata.obsm["protein_expression"].shape[1]
    n_latent = 10

    model = TOTALVI(adata, n_latent=n_latent)
    model.train(1, frequency=1, early_stopping_kwargs=None)
    assert len(model.history["elbo_train_set"]) == 2
    assert len(model.history["elbo_test_set"]) == 2
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
    adata2.uns["_scvi"]["categorical_mappings"]["_scvi_labels"]["mapping"] = np.array(
        ["undefined_1", "undefined_0", "undefined_2"]
    )
    with pytest.raises(ValueError):
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
    model.differential_expression(idx1=[0, 1, 2], idx2=[3, 4, 5])
    model.differential_expression(idx1=[0, 1, 2])
    model.differential_expression(groupby="labels")

    # test with missing proteins
    adata = scvi.data.pbmcs_10x_cite_seq(save_path=save_path, protein_join="outer")
    model = TOTALVI(adata)
    assert model.model.protein_batch_mask is not None
    model.train(2, train_size=0.5)


def test_scvi_online_update(save_path):
    n_latent = 5
    adata1 = synthetic_iid()
    model = SCVI(adata1, n_latent=n_latent)
    model.train(1, frequency=1)
    dir_path = os.path.join(save_path, "saved_model/")
    model.save(dir_path, overwrite=True)

    adata2 = synthetic_iid(run_setup_anndata=False)
    new_b = [2, 3]
    adata2.obs["batch"] = pd.Categorical(new_b[i] for i in adata2.obs.batch)

    model2 = SCVI.load_query_data(adata2, dir_path)
    model2.train(n_epochs=1)
    model2.get_latent_representation()

    # encoder linear layer equal
    one = model.model.z_encoder.encoder.fc_layers[0][0].weight.detach().numpy()
    two = model2.model.z_encoder.encoder.fc_layers[0][0].weight.detach().numpy()
    np.testing.assert_allclose(one, two, atol=1e-07)
    assert (
        np.sum(model2.model.z_encoder.encoder.fc_layers[0][0].weight.grad.numpy()) == 0
    )
    # dispersion
    assert model2.model.px_r.requires_grad is False
    # library encoder linear layer
    assert model2.model.l_encoder.encoder.fc_layers[0][0].weight.requires_grad is True
    # batch norm weight in encoder layer
    assert model2.model.z_encoder.encoder.fc_layers[0][1].weight.requires_grad is False
    # 5 for n_latent, 4 for batches
    assert model2.model.decoder.px_decoder.fc_layers[0][0].weight.shape[1] == 9

    # test options
    adata1 = synthetic_iid()
    model = SCVI(adata1, n_latent=n_latent, n_layers=2, encode_covariates=True)
    model.train(1, frequency=1)
    dir_path = os.path.join(save_path, "saved_model/")
    model.save(dir_path, overwrite=True)

    adata2 = synthetic_iid(run_setup_anndata=False)
    new_b = [2, 3]
    adata2.obs["batch"] = pd.Categorical(new_b[i] for i in adata2.obs.batch)

    model2 = SCVI.load_query_data(adata2, dir_path, freeze_expression=True)
    model2.train(n_epochs=1)
    model2.get_latent_representation()
    grad = model2.model.z_encoder.encoder.fc_layers[0][0].weight.grad.numpy()
    # expression part has zero grad
    assert np.sum(grad[:, :-4]) == 0
    # categorical part has non-zero grad
    assert np.sum(grad[:, -4:]) != 0

    # do not freeze expression
    model3 = SCVI.load_query_data(
        adata2, dir_path, freeze_expression=False, freeze_batchnorm=True
    )
    model3.train(n_epochs=1)
    model3.get_latent_representation()
    assert model3.model.z_encoder.encoder.fc_layers[0][1].momentum == 0
    grad = model3.model.z_encoder.encoder.fc_layers[0][0].weight.grad.numpy()
    # linear layer weight in encoder layer has non-zero grad
    assert np.sum(grad[:, :-4]) != 0

    # do not freeze batchnorm
    model3 = SCVI.load_query_data(adata2, dir_path, freeze_batchnorm=False)
    model3.train(n_epochs=1)
    model3.get_latent_representation()


def test_scanvi_online_update(save_path):
    n_latent = 5
    adata1 = synthetic_iid(run_setup_anndata=False)
    new_labels = adata1.obs.labels.to_numpy()
    new_labels[0] = "Unknown"
    adata1.obs["labels"] = pd.Categorical(new_labels)
    setup_anndata(adata1, batch_key="batch", labels_key="labels")
    model = SCANVI(adata1, "Unknown", n_latent=n_latent, encode_covariates=True)
    model.train(n_epochs_unsupervised=1, n_epochs_semisupervised=1, frequency=1)
    dir_path = os.path.join(save_path, "saved_model/")
    model.save(dir_path, overwrite=True)

    adata2 = synthetic_iid(run_setup_anndata=False)
    new_b = [2, 3]
    adata2.obs["batch"] = pd.Categorical(new_b[i] for i in adata2.obs.batch)
    adata2.obs["labels"] = "Unknown"

    model = SCANVI.load_query_data(adata2, dir_path, freeze_batchnorm=True)
    model.train(
        n_epochs_unsupervised=1, n_epochs_semisupervised=1, train_base_model=False
    )
    model.get_latent_representation()
    model.predict()


def test_totalvi_online_update(save_path):
    # basic case
    n_latent = 5
    adata1 = synthetic_iid()
    model = TOTALVI(adata1, n_latent=n_latent)
    model.train(1, frequency=1)
    dir_path = os.path.join(save_path, "saved_model/")
    model.save(dir_path, overwrite=True)

    adata2 = synthetic_iid(run_setup_anndata=False)
    new_b = [2, 3]
    adata2.obs["batch"] = pd.Categorical(new_b[i] for i in adata2.obs.batch)

    model = TOTALVI.load_query_data(adata2, dir_path)
    assert model.model.background_pro_alpha.requires_grad is True
    model.train(n_epochs=1)
    model.get_latent_representation()

    # batch 3 has no proteins
    adata2 = synthetic_iid(run_setup_anndata=False)
    new_b = [2, 3]
    adata2.obs["batch"] = pd.Categorical(new_b[i] for i in adata2.obs.batch)
    adata2.obsm["protein_expression"][adata2.obs.batch == 3] = 0

    model = TOTALVI.load_query_data(adata2, dir_path)
    model.model.protein_batch_mask[2]
    model.model.protein_batch_mask[3]
    model.train(n_epochs=1)
    model.get_latent_representation()
