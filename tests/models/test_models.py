import os
import tarfile

import anndata
import numpy as np
import pytest
from pytorch_lightning.callbacks import LearningRateMonitor
from scipy.sparse import csr_matrix
from torch.nn import Softplus

import scvi
from scvi.data import setup_anndata, synthetic_iid, transfer_anndata_setup
from scvi.data._built_in_data._download import _download
from scvi.dataloaders import (
    AnnDataLoader,
    DataSplitter,
    DeviceBackedDataSplitter,
    SemiSupervisedDataLoader,
    SemiSupervisedDataSplitter,
)
from scvi.model import (
    AUTOZI,
    PEAKVI,
    SCANVI,
    SCVI,
    TOTALVI,
    CondSCVI,
    DestVI,
    LinearSCVI,
)
from scvi.train import TrainingPlan, TrainRunner


def test_scvi(save_path):
    n_latent = 5
    adata = synthetic_iid()
    model = SCVI(adata, n_latent=n_latent)
    model.train(1, check_val_every_n_epoch=1, train_size=0.5)

    model = SCVI(adata, n_latent=n_latent, var_activation=Softplus())
    model.train(1, check_val_every_n_epoch=1, train_size=0.5)

    # tests __repr__
    print(model)

    assert model.is_trained is True
    z = model.get_latent_representation()
    assert z.shape == (adata.shape[0], n_latent)
    assert len(model.history["elbo_train"]) == 1
    model.get_elbo()
    model.get_marginal_ll(n_mc_samples=3)
    model.get_reconstruction_error()
    model.get_normalized_expression(transform_batch="batch_1")

    adata2 = synthetic_iid()
    model.get_elbo(adata2)
    model.get_marginal_ll(adata2, n_mc_samples=3)
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

    # test get_likelihood_parameters() when dispersion=='gene-cell'
    model = SCVI(adata, dispersion="gene-cell")
    model.get_likelihood_parameters()

    # test train callbacks work
    a = synthetic_iid()
    m = scvi.model.SCVI(a)
    lr_monitor = LearningRateMonitor()
    m.train(
        callbacks=[lr_monitor],
        max_epochs=10,
        log_every_n_steps=1,
        plan_kwargs={"reduce_lr_on_plateau": True},
    )
    assert "lr-Adam" in m.history.keys()


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
    model.get_marginal_ll(n_mc_samples=3)
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

    for cls in [SCVI, LinearSCVI, TOTALVI, PEAKVI]:
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
    a = scvi.data.synthetic_iid()
    m = scvi.model.SCVI.load(pretrained_scvi_path, a)
    m.train(1)
    pretrained_totalvi_path = os.path.join(save_path, "testing_models/080_totalvi")
    m = scvi.model.TOTALVI.load(pretrained_totalvi_path, a)
    m.train(1)


def test_backed_anndata_scvi(save_path):
    adata = scvi.data.synthetic_iid()
    path = os.path.join(save_path, "test_data.h5ad")
    adata.write_h5ad(path)
    adata = anndata.read_h5ad(path, backed="r+")
    setup_anndata(adata, batch_key="batch")

    model = SCVI(adata, n_latent=5)
    model.train(1, train_size=0.5)
    assert model.is_trained is True
    z = model.get_latent_representation()
    assert z.shape == (adata.shape[0], 5)
    model.get_elbo()


def test_ann_dataloader():
    a = scvi.data.synthetic_iid()

    # test that batch sampler drops the last batch if it has less than 3 cells
    assert a.n_obs == 400
    adl = AnnDataLoader(a, batch_size=397, drop_last=3)
    assert len(adl) == 2
    for i, x in enumerate(adl):
        pass
    assert i == 1
    adl = AnnDataLoader(a, batch_size=398, drop_last=3)
    assert len(adl) == 1
    for i, x in enumerate(adl):
        pass
    assert i == 0
    with pytest.raises(ValueError):
        AnnDataLoader(a, batch_size=1, drop_last=2)


def test_semisupervised_dataloader():
    # test label resampling
    n_samples_per_label = 10
    a = synthetic_iid()
    dl = SemiSupervisedDataLoader(
        a,
        indices=np.arange(a.n_obs),
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


def test_data_splitter():
    a = synthetic_iid()
    # test leaving validataion_size empty works
    ds = DataSplitter(a, train_size=0.4)
    ds.setup()
    # check the number of indices
    _, _, _ = ds.train_dataloader(), ds.val_dataloader(), ds.test_dataloader()
    n_train_idx = len(ds.train_idx)
    n_validation_idx = len(ds.val_idx) if ds.val_idx is not None else 0
    n_test_idx = len(ds.test_idx) if ds.test_idx is not None else 0

    assert n_train_idx + n_validation_idx + n_test_idx == a.n_obs
    assert np.isclose(n_train_idx / a.n_obs, 0.4)
    assert np.isclose(n_validation_idx / a.n_obs, 0.6)
    assert np.isclose(n_test_idx / a.n_obs, 0)

    # test test size
    ds = DataSplitter(a, train_size=0.4, validation_size=0.3)
    ds.setup()
    # check the number of indices
    _, _, _ = ds.train_dataloader(), ds.val_dataloader(), ds.test_dataloader()
    n_train_idx = len(ds.train_idx)
    n_validation_idx = len(ds.val_idx) if ds.val_idx is not None else 0
    n_test_idx = len(ds.test_idx) if ds.test_idx is not None else 0

    assert n_train_idx + n_validation_idx + n_test_idx == a.n_obs
    assert np.isclose(n_train_idx / a.n_obs, 0.4)
    assert np.isclose(n_validation_idx / a.n_obs, 0.3)
    assert np.isclose(n_test_idx / a.n_obs, 0.3)

    # test that 0 < train_size <= 1
    with pytest.raises(ValueError):
        ds = DataSplitter(a, train_size=2)
        ds.setup()
        ds.train_dataloader()
    with pytest.raises(ValueError):
        ds = DataSplitter(a, train_size=-2)
        ds.setup()
        ds.train_dataloader()

    # test that 0 <= validation_size < 1
    with pytest.raises(ValueError):
        ds = DataSplitter(a, train_size=0.1, validation_size=1)
        ds.setup()
        ds.val_dataloader()
    with pytest.raises(ValueError):
        ds = DataSplitter(a, train_size=0.1, validation_size=-1)
        ds.setup()
        ds.val_dataloader()

    # test that train_size + validation_size <= 1
    with pytest.raises(ValueError):
        ds = DataSplitter(a, train_size=1, validation_size=0.1)
        ds.setup()
        ds.train_dataloader()
        ds.val_dataloader()


def test_device_backed_data_splitter():
    a = synthetic_iid()
    # test leaving validataion_size empty works
    ds = DeviceBackedDataSplitter(a, train_size=1.0, use_gpu=None)
    ds.setup()
    train_dl = ds.train_dataloader()
    ds.val_dataloader()
    assert len(next(iter(train_dl))["X"]) == a.shape[0]

    model = SCVI(a, n_latent=5)
    training_plan = TrainingPlan(model.module, len(ds.train_idx))
    runner = TrainRunner(
        model,
        training_plan=training_plan,
        data_splitter=ds,
        max_epochs=1,
        use_gpu=None,
    )
    runner()


def test_semisupervised_data_splitter():
    a = synthetic_iid()
    ds = SemiSupervisedDataSplitter(a, "asdf")
    ds.setup()
    # check the number of indices
    _, _, _ = ds.train_dataloader(), ds.val_dataloader(), ds.test_dataloader()
    n_train_idx = len(ds.train_idx)
    n_validation_idx = len(ds.val_idx) if ds.val_idx is not None else 0
    n_test_idx = len(ds.test_idx) if ds.test_idx is not None else 0

    assert n_train_idx + n_validation_idx + n_test_idx == a.n_obs
    assert np.isclose(n_train_idx / a.n_obs, 0.9)
    assert np.isclose(n_validation_idx / a.n_obs, 0.1)
    assert np.isclose(n_test_idx / a.n_obs, 0)

    # test mix of labeled and unlabeled data
    unknown_label = "label_0"
    ds = SemiSupervisedDataSplitter(a, unknown_label)
    ds.setup()
    _, _, _ = ds.train_dataloader(), ds.val_dataloader(), ds.test_dataloader()

    # check the number of indices
    n_train_idx = len(ds.train_idx)
    n_validation_idx = len(ds.val_idx) if ds.val_idx is not None else 0
    n_test_idx = len(ds.test_idx) if ds.test_idx is not None else 0
    assert n_train_idx + n_validation_idx + n_test_idx == a.n_obs
    assert np.isclose(n_train_idx / a.n_obs, 0.9, rtol=0.05)
    assert np.isclose(n_validation_idx / a.n_obs, 0.1, rtol=0.05)
    assert np.isclose(n_test_idx / a.n_obs, 0, rtol=0.05)

    # check that training indices have proper mix of labeled and unlabeled data
    labelled_idx = np.where(a.obs["labels"] != unknown_label)[0]
    unlabelled_idx = np.where(a.obs["labels"] == unknown_label)[0]
    # labeled training idx
    labeled_train_idx = [i for i in ds.train_idx if i in labelled_idx]
    # unlabeled training idx
    unlabeled_train_idx = [i for i in ds.train_idx if i in unlabelled_idx]
    n_labeled_idx = len(labelled_idx)
    n_unlabeled_idx = len(unlabelled_idx)
    # labeled vs unlabeled ratio in adata
    adata_ratio = n_unlabeled_idx / n_labeled_idx
    # labeled vs unlabeled ratio in train set
    train_ratio = len(unlabeled_train_idx) / len(labeled_train_idx)
    assert np.isclose(adata_ratio, train_ratio, atol=0.05)


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
    assert "classification_loss_validation" in logged_keys
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

    # test mix of labeled and unlabeled data
    unknown_label = "label_0"
    a = scvi.data.synthetic_iid()
    scvi.data.setup_anndata(a, batch_key="batch", labels_key="labels")
    m = scvi.model.SCANVI(a, unknown_label)
    m.train(1, train_size=0.9)

    # test from_scvi_model
    a = scvi.data.synthetic_iid()
    m = scvi.model.SCVI(a, use_observed_lib_size=False)
    a2 = scvi.data.synthetic_iid()
    scanvi_model = scvi.model.SCANVI.from_scvi_model(m, "label_0", adata=a2)
    scanvi_model = scvi.model.SCANVI.from_scvi_model(
        m, "label_0", use_labels_groups=False
    )
    scanvi_model.train(1)


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


def test_autozi():
    data = synthetic_iid(n_batches=1)
    for disp_zi in ["gene", "gene-label"]:
        autozivae = AUTOZI(
            data,
            dispersion=disp_zi,
            zero_inflation=disp_zi,
        )
        autozivae.train(1, plan_kwargs=dict(lr=1e-2), check_val_every_n_epoch=1)
        assert len(autozivae.history["elbo_train"]) == 1
        assert len(autozivae.history["elbo_validation"]) == 1
        autozivae.get_elbo(indices=autozivae.validation_indices)
        autozivae.get_reconstruction_error(indices=autozivae.validation_indices)
        autozivae.get_marginal_ll(indices=autozivae.validation_indices, n_mc_samples=3)
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
    # model.get_likelihood_parameters()

    model.get_elbo(indices=model.validation_indices)
    model.get_marginal_ll(indices=model.validation_indices, n_mc_samples=3)
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
    assert model.module.protein_batch_mask is not None
    model.train(1, train_size=0.5)


def test_multiple_covariates(save_path):
    adata = synthetic_iid()
    adata.obs["cont1"] = np.random.normal(size=(adata.shape[0],))
    adata.obs["cont2"] = np.random.normal(size=(adata.shape[0],))
    adata.obs["cat1"] = np.random.randint(0, 5, size=(adata.shape[0],))
    adata.obs["cat2"] = np.random.randint(0, 5, size=(adata.shape[0],))
    setup_anndata(
        adata,
        batch_key="batch",
        labels_key="labels",
        protein_expression_obsm_key="protein_expression",
        protein_names_uns_key="protein_names",
        continuous_covariate_keys=["cont1", "cont2"],
        categorical_covariate_keys=["cat1", "cat2"],
    )

    m = SCVI(adata)
    m.train(1)

    m = SCANVI(adata, unlabeled_category="Unknown")
    m.train(1)

    m = TOTALVI(adata)
    m.train(1)


def test_peakvi():
    data = synthetic_iid()
    vae = PEAKVI(
        data,
        model_depth=False,
    )
    vae.train(1, save_best=False)
    vae = PEAKVI(
        data,
        region_factors=False,
    )
    vae.train(1, save_best=False)
    vae = PEAKVI(
        data,
    )
    vae.train(3)
    vae.get_elbo(indices=vae.validation_indices)
    vae.get_accessibility_estimates()
    vae.get_accessibility_estimates(normalize_cells=True)
    vae.get_accessibility_estimates(normalize_regions=True)
    vae.get_library_size_factors()
    vae.get_region_factors()
    vae.get_reconstruction_error(indices=vae.validation_indices)
    vae.get_latent_representation()
    vae.differential_accessibility(groupby="labels", group1="label_1")


def test_condscvi(save_path):
    dataset = synthetic_iid(n_labels=5)
    model = CondSCVI(dataset)
    model.train(1, train_size=1)
    model.get_latent_representation()
    model.get_vamp_prior(dataset)

    model = CondSCVI(dataset, weight_obs=True)
    model.train(1, train_size=1)
    model.get_latent_representation()
    model.get_vamp_prior(dataset)


def test_destvi(save_path):
    # Step1 learn CondSCVI
    n_latent = 2
    n_labels = 5
    n_layers = 2
    dataset = synthetic_iid(n_labels=n_labels)
    sc_model = CondSCVI(dataset, n_latent=n_latent, n_layers=n_layers)
    sc_model.train(1, train_size=1)

    # step 2 learn destVI with multiple amortization scheme

    for amor_scheme in ["both", "none", "proportion", "latent"]:
        spatial_model = DestVI.from_rna_model(
            dataset,
            sc_model,
            amortization=amor_scheme,
        )
        spatial_model.train(max_epochs=1)
        assert not np.isnan(spatial_model.history["elbo_train"].values[0][0])

        assert spatial_model.get_proportions().shape == (dataset.n_obs, n_labels)
        assert spatial_model.get_gamma(return_numpy=True).shape == (
            dataset.n_obs,
            n_latent,
            n_labels,
        )

        assert spatial_model.get_scale_for_ct("label_0", np.arange(50)).shape == (
            50,
            dataset.n_vars,
        )
