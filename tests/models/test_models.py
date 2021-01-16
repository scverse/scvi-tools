import os
import pandas as pd
import numpy as np
import pytest

import scvi
from scvi.data import synthetic_iid, transfer_anndata_setup, setup_anndata
from scvi.model import SCVI, SCANVI, GIMVI, TOTALVI, LinearSCVI, AUTOZI
from scvi.dataloaders import SemiSupervisedDataLoader

from scipy.sparse import csr_matrix


def test_scvi(save_path):
    n_latent = 5
    adata = synthetic_iid()
    model = SCVI(adata, n_latent=n_latent)
    model.train(1, check_val_every_n_epoch=1, train_size=0.5)

    # tests __repr__
    print(model)

    assert model.is_trained is True
    z = model.get_latent_representation()
    assert z.shape == (adata.shape[0], n_latent)
    assert len(model.history["elbo_train"]) == 1
    model.get_elbo()
    model.get_marginal_ll()
    model.get_reconstruction_error()
    model.get_normalized_expression(transform_batch="batch_1")

    adata2 = synthetic_iid()
    model.get_elbo(adata2)
    model.get_marginal_ll(adata2)
    model.get_reconstruction_error(adata2)
    latent = model.get_latent_representation(adata2, indices=[1, 2, 3])
    assert latent.shape == (3, n_latent)
    denoised = model.get_normalized_expression(adata2)
    assert denoised.shape == adata.shape

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
        ["label_4", "label_0", "label_2"]
    )
    with pytest.raises(ValueError):
        model.get_elbo(adata2)

    # test that same mapping different order doesn't raise error
    adata = synthetic_iid()
    adata2 = synthetic_iid(run_setup_anndata=False)
    transfer_anndata_setup(adata, adata2)
    adata2.uns["_scvi"]["categorical_mappings"]["_scvi_labels"]["mapping"] = np.array(
        ["label_1", "label_0", "label_2"]
    )
    model.get_elbo(adata2)  # should automatically transfer setup

    # test mismatched categories raises ValueError
    adata2 = synthetic_iid(run_setup_anndata=False)
    adata2.obs.labels.cat.rename_categories(["a", "b", "c"], inplace=True)
    with pytest.raises(ValueError):
        model.get_elbo(adata2)

    # test differential expression
    model.differential_expression(groupby="labels", group1="label_1")
    model.differential_expression(
        groupby="labels", group1="label_1", group2="label_2", mode="change"
    )
    model.differential_expression(groupby="labels")
    model.differential_expression(idx1=[0, 1, 2], idx2=[3, 4, 5])
    model.differential_expression(idx1=[0, 1, 2])

    # transform batch works with all different types
    a = synthetic_iid(run_setup_anndata=False)
    batch = np.zeros(a.n_obs)
    batch[:64] += 1
    a.obs["batch"] = batch
    setup_anndata(a, batch_key="batch")
    m = SCVI(a)
    m.train(1, train_size=0.5)
    m.get_normalized_expression(transform_batch=1)
    m.get_normalized_expression(transform_batch=[0, 1])


def test_scvi_sparse(save_path):
    n_latent = 5
    adata = synthetic_iid(run_setup_anndata=False)
    adata.X = csr_matrix(adata.X)
    setup_anndata(adata)
    model = SCVI(adata, n_latent=n_latent)
    model.train(1, train_size=0.5)
    assert model.is_trained is True
    z = model.get_latent_representation()
    assert z.shape == (adata.shape[0], n_latent)
    model.get_elbo()
    model.get_marginal_ll()
    model.get_reconstruction_error()
    model.get_normalized_expression()
    model.differential_expression(groupby="labels", group1="label_1")


def test_saving_and_loading(save_path):
    def test_save_load_model(cls, adata, save_path):
        model = cls(adata, latent_distribution="normal")
        model.train(1, train_size=0.2)
        z1 = model.get_latent_representation(adata)
        test_idx1 = model.validation_indices
        model.save(save_path, overwrite=True, save_anndata=True)
        model = cls.load(save_path)
        model.get_latent_representation()
        tmp_adata = scvi.data.synthetic_iid(n_genes=200)
        with pytest.raises(ValueError):
            cls.load(save_path, tmp_adata)
        model = cls.load(save_path, adata)
        z2 = model.get_latent_representation()
        test_idx2 = model.validation_indices
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
    model.train(1, train_size=0.5)
    ab1 = model.get_alphas_betas()
    model.save(save_path, overwrite=True, save_anndata=True)
    model = AUTOZI.load(save_path)
    model.get_latent_representation()
    tmp_adata = scvi.data.synthetic_iid(n_genes=200)
    with pytest.raises(ValueError):
        AUTOZI.load(save_path, tmp_adata)
    model = AUTOZI.load(save_path, adata)
    ab2 = model.get_alphas_betas()
    np.testing.assert_array_equal(ab1["alpha_posterior"], ab2["alpha_posterior"])
    np.testing.assert_array_equal(ab1["beta_posterior"], ab2["beta_posterior"])
    assert model.is_trained is True

    # SCANVI
    model = SCANVI(adata, "label_0")
    model.train(max_epochs=1, train_size=0.5)
    p1 = model.predict()
    model.save(save_path, overwrite=True, save_anndata=True)
    model = SCANVI.load(save_path)
    model.get_latent_representation()
    tmp_adata = scvi.data.synthetic_iid(n_genes=200)
    with pytest.raises(ValueError):
        SCANVI.load(save_path, tmp_adata)
    model = SCANVI.load(save_path, adata)
    p2 = model.predict()
    np.testing.assert_array_equal(p1, p2)
    assert model.is_trained is True

    # GIMVI
    model = GIMVI(adata, adata)
    model.train(1, train_size=0.5)
    z1 = model.get_latent_representation([adata])
    z2 = model.get_latent_representation([adata])
    np.testing.assert_array_equal(z1, z2)
    model.save(save_path, overwrite=True, save_anndata=True)
    model = GIMVI.load(save_path)
    model.get_latent_representation()
    tmp_adata = scvi.data.synthetic_iid(n_genes=200)
    with pytest.raises(ValueError):
        GIMVI.load(save_path, tmp_adata, tmp_adata)
    model = GIMVI.load(save_path, adata, adata)
    z2 = model.get_latent_representation([adata])
    np.testing.assert_array_equal(z1, z2)
    assert model.is_trained is True


def test_SemiSupervisedDataLoader():
    # test label resampling
    n_samples_per_label = 10
    a = synthetic_iid()
    dl = SemiSupervisedDataLoader(
        a,
        indices=np.arange(a.n_obs),
        labels_obs_key="labels",
        unlabeled_category="label_0",
        n_samples_per_label=n_samples_per_label,
    )
    labeled_dl_idx = dl.dataloaders[1].indices
    n_labels = 2
    assert len(labeled_dl_idx) == n_samples_per_label * n_labels
    dl.resample_labels()
    resampled_labeled_dl_idx = dl.dataloaders[1].indices
    assert len(resampled_labeled_dl_idx) == n_samples_per_label * n_labels
    # check labeled indices was actually resampled
    assert np.sum(labeled_dl_idx == resampled_labeled_dl_idx) != len(labeled_dl_idx)


def test_scanvi(save_path):
    adata = synthetic_iid()
    model = SCANVI(adata, "label_0", n_latent=10)
    model.train(1, train_size=0.5, check_val_every_n_epoch=1)
    logged_keys = model.history.keys()
    assert "elbo_validation" in logged_keys
    assert "reconstruction_loss_validation" in logged_keys
    assert "kl_local_validation" in logged_keys
    assert "elbo_train" in logged_keys
    assert "reconstruction_loss_train" in logged_keys
    assert "kl_local_train" in logged_keys
    adata2 = synthetic_iid()
    predictions = model.predict(adata2, indices=[1, 2, 3])
    assert len(predictions) == 3
    model.predict()
    model.predict(adata2, soft=True)
    model.predict(adata2, soft=True, indices=[1, 2, 3])
    model.get_normalized_expression(adata2)
    model.differential_expression(groupby="labels", group1="label_1")
    model.differential_expression(groupby="labels", group1="label_1", group2="label_2")

    # test that all data labeled runs
    unknown_label = "asdf"
    a = scvi.data.synthetic_iid()
    scvi.data.setup_anndata(a, batch_key="batch", labels_key="labels")
    m = scvi.model.SCANVI(a, unknown_label)
    m.train(1)

    # check the number of indices
    n_train_idx = len(m.train_indices)
    n_validation_idx = len(m.validation_indices)
    n_test_idx = len(m.test_indices)
    assert n_train_idx + n_validation_idx + n_test_idx == a.n_obs
    assert np.isclose(n_train_idx / a.n_obs, 0.9)
    assert np.isclose(n_validation_idx / a.n_obs, 0.1)
    assert np.isclose(n_test_idx / a.n_obs, 0)

    # test mix of labeled and unlabeled data
    unknown_label = "label_0"
    a = scvi.data.synthetic_iid()
    scvi.data.setup_anndata(a, batch_key="batch", labels_key="labels")
    m = scvi.model.SCANVI(a, unknown_label)
    m.train(1, train_size=0.9)
    # check the number of indices
    n_train_idx = len(m.train_indices)
    n_validation_idx = len(m.validation_indices)
    n_test_idx = len(m.test_indices)
    assert n_train_idx + n_validation_idx + n_test_idx == a.n_obs
    assert np.isclose(n_train_idx / a.n_obs, 0.9, rtol=0.05)
    assert np.isclose(n_validation_idx / a.n_obs, 0.1, rtol=0.05)
    assert np.isclose(n_test_idx / a.n_obs, 0, rtol=0.05)

    # check that training indices have proper mix of labeled and unlabeled data
    labelled_idx = np.where(a.obs["labels"] != unknown_label)[0]
    unlabelled_idx = np.where(a.obs["labels"] == unknown_label)[0]
    # labeled training idx
    labeled_train_idx = [i for i in m.train_indices if i in labelled_idx]
    # unlabeled training idx
    unlabeled_train_idx = [i for i in m.train_indices if i in unlabelled_idx]
    n_labeled_idx = len(m._labeled_indices)
    n_unlabeled_idx = len(m._unlabeled_indices)
    # labeled vs unlabeled ratio in adata
    adata_ratio = n_unlabeled_idx / n_labeled_idx
    # labeled vs unlabeled ratio in train set
    train_ratio = len(unlabeled_train_idx) / len(labeled_train_idx)
    assert np.isclose(adata_ratio, train_ratio, atol=0.05)


def test_linear_scvi(save_path):
    adata = synthetic_iid()
    adata = adata[:, :10].copy()
    setup_anndata(adata)
    model = LinearSCVI(adata, n_latent=10)
    model.train(1, check_val_every_n_epoch=1, train_size=0.5)
    assert len(model.history["elbo_train"]) == 1
    assert len(model.history["elbo_validation"]) == 1
    model.get_loadings()
    model.differential_expression(groupby="labels", group1="label_1")
    model.differential_expression(groupby="labels", group1="label_1", group2="label_2")


def test_gimvi():
    adata_seq = synthetic_iid()
    adata_spatial = synthetic_iid()
    model = GIMVI(adata_seq, adata_spatial, n_latent=10)
    model.get_latent_representation()
    model.get_imputed_values()
    model.train(
        1, check_val_every_n_epoch=1, early_stopping_kwargs=None, train_size=0.5
    )

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
        autozivae.train(1, lr=1e-2, check_val_every_n_epoch=1)
        assert len(autozivae.history["elbo_train"]) == 1
        assert len(autozivae.history["elbo_test"]) == 1
        autozivae.get_elbo(indices=autozivae.validation_indices)
        autozivae.get_reconstruction_error(indices=autozivae.validation_indices)
        autozivae.get_marginal_ll(indices=autozivae.validation_indices)
        autozivae.get_alphas_betas()


def test_totalvi(save_path):
    adata = synthetic_iid()
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
    model.get_marginal_ll()
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
    # model.get_likelihood_parameters()

    model.get_elbo(indices=model.validation_indices)
    model.get_marginal_ll(indices=model.validation_indices)
    model.get_reconstruction_error(indices=model.validation_indices)

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
        ["label_1", "label_0", "label_8"]
    )
    with pytest.raises(ValueError):
        model.get_elbo(adata2)

    # test that same mapping different order is okay
    adata = synthetic_iid()
    adata2 = synthetic_iid(run_setup_anndata=False)
    transfer_anndata_setup(adata, adata2)
    adata2.uns["_scvi"]["categorical_mappings"]["_scvi_labels"]["mapping"] = np.array(
        ["label_1", "label_0", "label_2"]
    )
    model.get_elbo(adata2)  # should automatically transfer setup

    # test that we catch missing proteins
    adata2 = synthetic_iid(run_setup_anndata=False)
    del adata2.obsm["protein_expression"]
    with pytest.raises(KeyError):
        model.get_elbo(adata2)
    model.differential_expression(groupby="labels", group1="label_1")
    model.differential_expression(groupby="labels", group1="label_1", group2="label_2")
    model.differential_expression(idx1=[0, 1, 2], idx2=[3, 4, 5])
    model.differential_expression(idx1=[0, 1, 2])
    model.differential_expression(groupby="labels")

    # test with missing proteins
    adata = scvi.data.pbmcs_10x_cite_seq(save_path=save_path, protein_join="outer")
    model = TOTALVI(adata)
    assert model.model.protein_batch_mask is not None
    model.train(1, train_size=0.5)


def single_pass_for_online_update(model):
    dl = model._make_scvi_dl(model.adata, indices=range(0, 10))
    for i_batch, tensors in enumerate(dl):
        _, _, scvi_loss = model.model(tensors)
    scvi_loss.loss.backward()


def test_scvi_online_update(save_path):
    n_latent = 5
    adata1 = synthetic_iid()
    model = SCVI(adata1, n_latent=n_latent)
    model.train(1, check_val_every_n_epoch=1)
    dir_path = os.path.join(save_path, "saved_model/")
    model.save(dir_path, overwrite=True)

    # also test subset var option
    adata2 = synthetic_iid(run_setup_anndata=False, n_genes=110)
    adata2.obs["batch"] = adata2.obs.batch.cat.rename_categories(["batch_2", "batch_3"])

    model2 = SCVI.load_query_data(adata2, dir_path, inplace_subset_query_vars=True)
    model2.train(max_epochs=1, vae_task_kwargs=dict(weight_decay=0.0))
    model2.get_latent_representation()

    # encoder linear layer equal
    one = (
        model.model.z_encoder.encoder.fc_layers[0][0]
        .weight.detach()
        .numpy()[:, : adata1.shape[1]]
    )
    two = (
        model2.model.z_encoder.encoder.fc_layers[0][0]
        .weight.detach()
        .numpy()[:, : adata1.shape[1]]
    )
    np.testing.assert_equal(one, two)
    assert (
        np.sum(
            model2.model.z_encoder.encoder.fc_layers[0][0].weight.grad.numpy()[
                :, : adata1.shape[1]
            ]
        )
        == 0
    )
    # dispersion
    assert model2.model.px_r.requires_grad is False
    # library encoder linear layer
    assert model2.model.l_encoder.encoder.fc_layers[0][0].weight.requires_grad is True
    # 5 for n_latent, 4 for batches
    assert model2.model.decoder.px_decoder.fc_layers[0][0].weight.shape[1] == 9

    # test options
    adata1 = synthetic_iid()
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

    adata2 = synthetic_iid(run_setup_anndata=False)
    adata2.obs["batch"] = adata2.obs.batch.cat.rename_categories(["batch_2", "batch_3"])

    model2 = SCVI.load_query_data(adata2, dir_path, freeze_expression=True)
    model2.train(max_epochs=1, vae_task_kwargs=dict(weight_decay=0.0))
    # deactivate no grad decorator
    model2.get_latent_representation()
    # pytorch lightning zeros the grad, so this will get a grad to inspect
    single_pass_for_online_update(model2)
    grad = model2.model.z_encoder.encoder.fc_layers[0][0].weight.grad.numpy()
    # expression part has zero grad
    assert np.sum(grad[:, :-4]) == 0
    # categorical part has non-zero grad
    assert np.sum(grad[:, -4:]) != 0

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
    assert model3.model.z_encoder.encoder.fc_layers[0][1].momentum == 0
    # batch norm weight in encoder layer
    assert model3.model.z_encoder.encoder.fc_layers[0][1].weight.requires_grad is False
    single_pass_for_online_update(model3)
    grad = model3.model.z_encoder.encoder.fc_layers[0][0].weight.grad.numpy()
    # linear layer weight in encoder layer has non-zero grad
    assert np.sum(grad[:, :-4]) != 0
    grad = model3.model.decoder.px_decoder.fc_layers[0][0].weight.grad.numpy()
    # linear layer weight in decoder layer has non-zero grad
    assert np.sum(grad[:, :-4]) != 0

    # do not freeze batchnorm
    model3 = SCVI.load_query_data(adata2, dir_path, freeze_batchnorm_encoder=False)
    model3.train(max_epochs=1)
    model3.get_latent_representation()


def test_scanvi_online_update(save_path):
    # ref has semi-observed labels
    n_latent = 5
    adata1 = synthetic_iid(run_setup_anndata=False)
    new_labels = adata1.obs.labels.to_numpy()
    new_labels[0] = "Unknown"
    adata1.obs["labels"] = pd.Categorical(new_labels)
    setup_anndata(adata1, batch_key="batch", labels_key="labels")
    model = SCANVI(adata1, "Unknown", n_latent=n_latent, encode_covariates=True)
    model.train(max_epochs=1, check_val_every_n_epoch=1)
    dir_path = os.path.join(save_path, "saved_model/")
    model.save(dir_path, overwrite=True)

    adata2 = synthetic_iid(run_setup_anndata=False)
    adata2.obs["batch"] = adata2.obs.batch.cat.rename_categories(["batch_2", "batch_3"])
    adata2.obs["labels"] = "Unknown"

    model = SCANVI.load_query_data(adata2, dir_path, freeze_batchnorm_encoder=True)
    model.train(max_epochs=1)
    model.get_latent_representation()
    model.predict()

    # ref has fully-observed labels
    n_latent = 5
    adata1 = synthetic_iid(run_setup_anndata=False)
    new_labels = adata1.obs.labels.to_numpy()
    adata1.obs["labels"] = pd.Categorical(new_labels)
    setup_anndata(adata1, batch_key="batch", labels_key="labels")
    model = SCANVI(adata1, "Unknown", n_latent=n_latent, encode_covariates=True)
    model.train(max_epochs=1, check_val_every_n_epoch=1)
    dir_path = os.path.join(save_path, "saved_model/")
    model.save(dir_path, overwrite=True)

    # query has one new label
    adata2 = synthetic_iid(run_setup_anndata=False)
    adata2.obs["batch"] = adata2.obs.batch.cat.rename_categories(["batch_2", "batch_3"])
    new_labels = adata2.obs.labels.to_numpy()
    new_labels[0] = "Unknown"
    adata2.obs["labels"] = pd.Categorical(new_labels)

    model2 = SCANVI.load_query_data(adata2, dir_path, freeze_batchnorm_encoder=True)
    model2._unlabeled_indices = np.arange(adata2.n_obs)
    model2._labeled_indices = []
    model2.train(max_epochs=1, task_kwargs=dict(weight_decay=0.0))
    model2.get_latent_representation()
    model2.predict()

    # test classifier frozen
    class_query_weight = (
        model2.model.classifier.classifier[0].fc_layers[0][0].weight.detach().numpy()
    )
    class_ref_weight = (
        model.model.classifier.classifier[0].fc_layers[0][0].weight.detach().numpy()
    )
    # weight decay makes difference
    np.testing.assert_allclose(class_query_weight, class_ref_weight, atol=1e-07)

    # test classifier unfrozen
    model2 = SCANVI.load_query_data(adata2, dir_path, freeze_classifier=False)
    model2._unlabeled_indices = np.arange(adata2.n_obs)
    model2._labeled_indices = []
    model2.train(max_epochs=1)
    class_query_weight = (
        model2.model.classifier.classifier[0].fc_layers[0][0].weight.detach().numpy()
    )
    class_ref_weight = (
        model.model.classifier.classifier[0].fc_layers[0][0].weight.detach().numpy()
    )
    with pytest.raises(AssertionError):
        np.testing.assert_allclose(class_query_weight, class_ref_weight, atol=1e-07)

    # test saving and loading of online scanvi
    a = scvi.data.synthetic_iid(run_setup_anndata=False)
    ref = a[a.obs["labels"] != "label_2"].copy()  # only has labels 0 and 1
    scvi.data.setup_anndata(ref, batch_key="batch", labels_key="labels")
    m = scvi.model.SCANVI(ref, "label_2")
    m.train(1, 1)
    m.save(save_path, overwrite=True)
    query = a[a.obs["labels"] != "label_0"].copy()
    query = scvi.data.synthetic_iid()  # has labels 0 and 2. 2 is unknown
    m_q = scvi.model.SCANVI.load_query_data(query, save_path)
    m_q.save(save_path, overwrite=True)
    m_q = scvi.model.SCANVI.load(save_path, query)
    m_q.predict()
    m_q.get_elbo()


def test_totalvi_online_update(save_path):
    # basic case
    n_latent = 5
    adata1 = synthetic_iid()
    model = TOTALVI(adata1, n_latent=n_latent, use_batch_norm="decoder")
    model.train(1, check_val_every_n_epoch=1)
    dir_path = os.path.join(save_path, "saved_model/")
    model.save(dir_path, overwrite=True)

    adata2 = synthetic_iid(run_setup_anndata=False)
    adata2.obs["batch"] = adata2.obs.batch.cat.rename_categories(["batch_2", "batch_3"])

    model2 = TOTALVI.load_query_data(adata2, dir_path)
    assert model2.model.background_pro_alpha.requires_grad is True
    model2.train(max_epochs=1)
    model2.get_latent_representation()

    # batch 3 has no proteins
    adata2 = synthetic_iid(run_setup_anndata=False)
    adata2.obs["batch"] = adata2.obs.batch.cat.rename_categories(["batch_2", "batch_3"])
    adata2.obsm["protein_expression"][adata2.obs.batch == "batch_3"] = 0

    # load from model in memory
    model3 = TOTALVI.load_query_data(adata2, model)
    model3.model.protein_batch_mask[2]
    model3.model.protein_batch_mask[3]
    model3.train(max_epochs=1)
    model3.get_latent_representation()
