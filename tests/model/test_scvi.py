import inspect
import os
import tarfile
from unittest import mock

import anndata
import numpy as np
import pandas as pd
import pytest
from lightning.pytorch.callbacks import LearningRateMonitor
from scipy.sparse import csr_matrix
from torch.nn import Softplus

import scvi
from scvi.data import _constants, synthetic_iid
from scvi.data._compat import LEGACY_REGISTRY_KEY_MAP, registry_from_setup_dict
from scvi.data._download import _download
from scvi.dataloaders import (
    DeviceBackedDataSplitter,
)
from scvi.model import (
    SCANVI,
    SCVI,
    TOTALVI,
)
from scvi.model.utils import mde
from scvi.train import TrainingPlan, TrainRunner

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


def test_scvi(save_path):
    n_latent = 5

    # Test with size factor.
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
    m = scvi.model.SCVI(a)
    lr_monitor = LearningRateMonitor()
    m.train(
        callbacks=[lr_monitor],
        max_epochs=10,
        check_val_every_n_epoch=1,
        log_every_n_steps=1,
        plan_kwargs={"reduce_lr_on_plateau": True},
    )
    assert "lr-Adam" in m.history.keys()


def test_scvi_get_latent_rep_backwards_compat():
    n_latent = 5

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


def test_scvi_get_feature_corr_backwards_compat():
    n_latent = 5

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


def test_scvi_sparse(save_path):
    n_latent = 5
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


def test_setting_adata_attr():
    n_latent = 5
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


def test_backed_anndata_scvi(save_path):
    adata = scvi.data.synthetic_iid()
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


def test_scanvi(save_path):
    adata = synthetic_iid()
    SCANVI.setup_anndata(
        adata,
        "labels",
        "label_0",
        batch_key="batch",
    )
    model = SCANVI(adata, n_latent=10)
    assert len(model._labeled_indices) == sum(adata.obs["labels"] != "label_0")
    assert len(model._unlabeled_indices) == sum(adata.obs["labels"] == "label_0")
    model.train(1, train_size=0.5, check_val_every_n_epoch=1)
    logged_keys = model.history.keys()
    assert "elbo_validation" in logged_keys
    assert "reconstruction_loss_validation" in logged_keys
    assert "kl_local_validation" in logged_keys
    assert "elbo_train" in logged_keys
    assert "reconstruction_loss_train" in logged_keys
    assert "kl_local_train" in logged_keys
    assert "validation_classification_loss" in logged_keys
    adata2 = synthetic_iid()
    predictions = model.predict(adata2, indices=[1, 2, 3])
    assert len(predictions) == 3
    model.predict()
    df = model.predict(adata2, soft=True)
    assert isinstance(df, pd.DataFrame)
    model.predict(adata2, soft=True, indices=[1, 2, 3])
    model.get_normalized_expression(adata2)
    model.differential_expression(groupby="labels", group1="label_1")
    model.differential_expression(groupby="labels", group1="label_1", group2="label_2")

    # test that all data labeled runs
    unknown_label = "asdf"
    a = scvi.data.synthetic_iid()
    scvi.model.SCANVI.setup_anndata(
        a,
        "labels",
        unknown_label,
        batch_key="batch",
    )
    m = scvi.model.SCANVI(a)
    m.train(1)

    # test mix of labeled and unlabeled data
    unknown_label = "label_0"
    a = scvi.data.synthetic_iid()
    scvi.model.SCANVI.setup_anndata(
        a,
        "labels",
        unknown_label,
        batch_key="batch",
    )
    m = scvi.model.SCANVI(a)
    m.train(1, train_size=0.9)

    # test from_scvi_model
    a = scvi.data.synthetic_iid()
    SCVI.setup_anndata(
        a,
        batch_key="batch",
    )
    m = SCVI(a, use_observed_lib_size=False)
    a2 = scvi.data.synthetic_iid()
    scanvi_model = scvi.model.SCANVI.from_scvi_model(
        m, "label_0", labels_key="labels", adata=a2
    )
    with pytest.raises(ValueError):
        scanvi_model = scvi.model.SCANVI.from_scvi_model(
            m, "label_0", labels_key=None, adata=a2
        )

    # make sure the state_dicts are different objects for the two models
    assert scanvi_model.module.state_dict() is not m.module.state_dict()
    scanvi_pxr = scanvi_model.module.state_dict().get("px_r", None)
    scvi_pxr = m.module.state_dict().get("px_r", None)
    assert scanvi_pxr is not None and scvi_pxr is not None
    assert scanvi_pxr is not scvi_pxr
    scanvi_model.train(1)

    # Test without label groups
    scanvi_model = scvi.model.SCANVI.from_scvi_model(
        m, "label_0", labels_key="labels", use_labels_groups=False
    )
    scanvi_model.train(1)

    # test from_scvi_model with size_factor
    a = scvi.data.synthetic_iid()
    a.obs["size_factor"] = np.random.randint(1, 5, size=(a.shape[0],))
    SCVI.setup_anndata(
        a, batch_key="batch", labels_key="labels", size_factor_key="size_factor"
    )
    m = SCVI(a, use_observed_lib_size=False)
    a2 = scvi.data.synthetic_iid()
    a2.obs["size_factor"] = np.random.randint(1, 5, size=(a2.shape[0],))
    scanvi_model = scvi.model.SCANVI.from_scvi_model(m, "label_0", adata=a2)
    scanvi_model.train(1)


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
