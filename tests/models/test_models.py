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
from scvi.dataloaders import (
    AnnDataLoader,
    DataSplitter,
    DeviceBackedDataSplitter,
    SemiSupervisedDataLoader,
    SemiSupervisedDataSplitter,
)
from scvi.model import (
    AUTOZI,
    MULTIVI,
    PEAKVI,
    SCANVI,
    SCVI,
    TOTALVI,
    CondSCVI,
    DestVI,
    LinearSCVI,
)
from scvi.train import TrainingPlan, TrainRunner
from scvi.utils import attrdict
from tests.dataset.utils import generic_setup_adata_manager, scanvi_setup_adata_manager

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


@pytest.mark.parametrize(
    "data",
    [
        scvi.data.synthetic_iid(200),
        scvi.data.synthetic_iid(200, sparse_format="csr_array"),
    ],
)
def test_ann_dataloader(data):
    adata_manager = generic_setup_adata_manager(
        data, batch_key="batch", labels_key="labels"
    )

    # test that batch sampler drops the last batch if it has less than 3 cells
    assert data.n_obs == 400
    adl = AnnDataLoader(adata_manager, batch_size=397, drop_last=True)
    assert len(adl) == 1
    for _i, _ in enumerate(adl):
        pass
    assert _i == 0
    adl = AnnDataLoader(adata_manager, batch_size=397, drop_last=False)
    assert len(adl) == 2
    for _i, _ in enumerate(adl):
        pass
    assert _i == 1
    adl = AnnDataLoader(adata_manager, batch_size=399, drop_last=False)
    assert len(adl) == 2
    for _i, loaded_data in enumerate(adl):
        _ = loaded_data
    assert _i == 1
    assert loaded_data["X"].shape[0] == 1


def test_semisupervised_dataloader():
    # test label resampling
    n_samples_per_label = 10
    a = synthetic_iid()
    adata_manager = scanvi_setup_adata_manager(
        a, labels_key="labels", unlabeled_category="label_0", batch_key="batch"
    )
    dl = SemiSupervisedDataLoader(
        adata_manager,
        indices=np.arange(a.n_obs),
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
    adata_manager = generic_setup_adata_manager(
        a, batch_key="batch", labels_key="labels"
    )
    # test leaving validataion_size empty works
    ds = DataSplitter(adata_manager, train_size=0.4)
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
    ds = DataSplitter(adata_manager, train_size=0.4, validation_size=0.3)
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
        ds = DataSplitter(adata_manager, train_size=2)
        ds.setup()
        ds.train_dataloader()
    with pytest.raises(ValueError):
        ds = DataSplitter(adata_manager, train_size=-2)
        ds.setup()
        ds.train_dataloader()

    # test that 0 <= validation_size < 1
    with pytest.raises(ValueError):
        ds = DataSplitter(adata_manager, train_size=0.1, validation_size=1)
        ds.setup()
        ds.val_dataloader()
    with pytest.raises(ValueError):
        ds = DataSplitter(adata_manager, train_size=0.1, validation_size=-1)
        ds.setup()
        ds.val_dataloader()

    # test that train_size + validation_size <= 1
    with pytest.raises(ValueError):
        ds = DataSplitter(adata_manager, train_size=1, validation_size=0.1)
        ds.setup()
        ds.train_dataloader()
        ds.val_dataloader()


def test_device_backed_data_splitter():
    a = synthetic_iid()
    SCVI.setup_anndata(a, batch_key="batch", labels_key="labels")
    model = SCVI(a, n_latent=5)
    adata_manager = model.adata_manager
    # test leaving validataion_size empty works
    ds = DeviceBackedDataSplitter(adata_manager, train_size=1.0)
    ds.setup()
    train_dl = ds.train_dataloader()
    ds.val_dataloader()
    loaded_x = next(iter(train_dl))["X"]
    assert len(loaded_x) == a.shape[0]
    np.testing.assert_array_equal(loaded_x.cpu().numpy(), a.X)

    training_plan = TrainingPlan(model.module)
    runner = TrainRunner(
        model,
        training_plan=training_plan,
        data_splitter=ds,
        max_epochs=1,
    )
    runner()


def test_semisupervised_data_splitter():
    a = synthetic_iid()
    adata_manager = scanvi_setup_adata_manager(
        a, labels_key="labels", unlabeled_category="asdf", batch_key="batch"
    )
    ds = SemiSupervisedDataSplitter(adata_manager)
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
    ds = SemiSupervisedDataSplitter(adata_manager)
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


def test_linear_scvi(save_path):
    adata = synthetic_iid()
    adata = adata[:, :10].copy()
    LinearSCVI.setup_anndata(adata)
    model = LinearSCVI(adata, n_latent=10)
    model.train(1, check_val_every_n_epoch=1, train_size=0.5)
    assert len(model.history["elbo_train"]) == 1
    assert len(model.history["elbo_validation"]) == 1
    model.get_loadings()
    model.differential_expression(groupby="labels", group1="label_1")
    model.differential_expression(groupby="labels", group1="label_1", group2="label_2")


def test_autozi():
    data = synthetic_iid(
        n_batches=1,
    )
    AUTOZI.setup_anndata(
        data,
        batch_key="batch",
        labels_key="labels",
    )

    for disp_zi in ["gene", "gene-label"]:
        autozivae = AUTOZI(
            data,
            dispersion=disp_zi,
            zero_inflation=disp_zi,
        )
        autozivae.train(1, plan_kwargs={"lr": 1e-2}, check_val_every_n_epoch=1)
        assert len(autozivae.history["elbo_train"]) == 1
        assert len(autozivae.history["elbo_validation"]) == 1
        autozivae.get_elbo(indices=autozivae.validation_indices)
        autozivae.get_reconstruction_error(indices=autozivae.validation_indices)
        autozivae.get_marginal_ll(indices=autozivae.validation_indices, n_mc_samples=3)
        autozivae.get_alphas_betas()

    # Model library size.
    for disp_zi in ["gene", "gene-label"]:
        autozivae = AUTOZI(
            data,
            dispersion=disp_zi,
            zero_inflation=disp_zi,
            use_observed_lib_size=False,
        )
        autozivae.train(1, plan_kwargs={"lr": 1e-2}, check_val_every_n_epoch=1)
        assert hasattr(autozivae.module, "library_log_means") and hasattr(
            autozivae.module, "library_log_vars"
        )
        assert len(autozivae.history["elbo_train"]) == 1
        assert len(autozivae.history["elbo_validation"]) == 1
        autozivae.get_elbo(indices=autozivae.validation_indices)
        autozivae.get_reconstruction_error(indices=autozivae.validation_indices)
        autozivae.get_marginal_ll(indices=autozivae.validation_indices, n_mc_samples=3)
        autozivae.get_alphas_betas()


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


def test_peakvi():
    data = synthetic_iid()
    PEAKVI.setup_anndata(
        data,
        batch_key="batch",
    )
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
    dataset = synthetic_iid(
        n_labels=5,
    )
    CondSCVI.setup_anndata(
        dataset,
        "labels",
    )
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
    dataset.obs["overclustering_vamp"] = list(range(dataset.n_obs))
    CondSCVI.setup_anndata(dataset, labels_key="labels")
    sc_model = CondSCVI(dataset, n_latent=n_latent, n_layers=n_layers)
    sc_model.train(1, train_size=1)

    # step 2 Check model setup
    DestVI.setup_anndata(dataset, layer=None)

    # Test clustering outside of get_vamp_prior

    # vamp_prior_p>n_largest_cluster to be successful.
    _ = DestVI.from_rna_model(dataset, sc_model, vamp_prior_p=dataset.n_obs)
    # vamp_prior_p<n_largest_cluster leads to value error.
    with pytest.raises(ValueError):
        _ = DestVI.from_rna_model(dataset, sc_model, vamp_prior_p=1)

    del dataset.obs["overclustering_vamp"]

    # step 3 learn destVI with multiple amortization scheme

    for amor_scheme in ["both", "none", "proportion", "latent"]:
        DestVI.setup_anndata(dataset, layer=None)
        # add l1_regularization to cell type proportions
        spatial_model = DestVI.from_rna_model(
            dataset, sc_model, amortization=amor_scheme, l1_reg=50
        )
        spatial_model.view_anndata_setup()
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


def test_multivi():
    data = synthetic_iid()
    MULTIVI.setup_anndata(
        data,
        batch_key="batch",
    )
    vae = MULTIVI(
        data,
        n_genes=50,
        n_regions=50,
    )
    vae.train(1, save_best=False)
    vae.train(1, adversarial_mixing=False)
    vae.train(3)
    vae.get_elbo(indices=vae.validation_indices)
    vae.get_accessibility_estimates()
    vae.get_accessibility_estimates(normalize_cells=True)
    vae.get_accessibility_estimates(normalize_regions=True)
    vae.get_normalized_expression()
    vae.get_library_size_factors()
    vae.get_region_factors()
    vae.get_reconstruction_error(indices=vae.validation_indices)
    vae.get_latent_representation()
    vae.differential_accessibility(groupby="labels", group1="label_1")
    vae.differential_expression(groupby="labels", group1="label_1")

    # Test with size factor
    data = synthetic_iid()
    data.obs["size_factor"] = np.random.randint(1, 5, size=(data.shape[0],))
    MULTIVI.setup_anndata(data, batch_key="batch", size_factor_key="size_factor")
    vae = MULTIVI(
        data,
        n_genes=50,
        n_regions=50,
    )
    vae.train(3)

    # Test with modality weights and penalties
    data = synthetic_iid()
    MULTIVI.setup_anndata(data, batch_key="batch")
    vae = MULTIVI(data, n_genes=50, n_regions=50, modality_weights="cell")
    vae.train(3)
    vae = MULTIVI(data, n_genes=50, n_regions=50, modality_weights="universal")
    vae.train(3)
    vae = MULTIVI(data, n_genes=50, n_regions=50, modality_penalty="MMD")
    vae.train(3)

    # Test with non-zero protein data
    data = synthetic_iid()
    MULTIVI.setup_anndata(
        data,
        batch_key="batch",
        protein_expression_obsm_key="protein_expression",
        protein_names_uns_key="protein_names",
    )
    vae = MULTIVI(
        data,
        n_genes=50,
        n_regions=50,
        modality_weights="cell",
    )
    assert vae.n_proteins == data.obsm["protein_expression"].shape[1]
    vae.train(3)
