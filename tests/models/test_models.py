import os
import pickle
import tarfile

import anndata
import numpy as np
import pandas as pd
import pytest
import torch
from pytorch_lightning.callbacks import LearningRateMonitor
from scipy.sparse import csr_matrix
from torch.nn import Softplus

import scvi
from scvi.data import _constants, synthetic_iid
from scvi.data._built_in_data._download import _download
from scvi.data._compat import LEGACY_REGISTRY_KEY_MAP, manager_from_setup_dict
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
    JaxSCVI,
    LinearSCVI,
)
from scvi.model.utils import mde
from scvi.train import TrainingPlan, TrainRunner
from tests.dataset.utils import generic_setup_adata_manager

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


def test_jax_scvi():
    n_latent = 5

    # Test with size factor.
    adata = synthetic_iid()
    JaxSCVI.setup_anndata(
        adata,
        batch_key="batch",
    )
    model = JaxSCVI(adata, n_latent=n_latent)
    model.train(2, train_size=0.5, check_val_every_n_epoch=1)
    model.get_latent_representation()

    model = JaxSCVI(adata, n_latent=n_latent, gene_likelihood="poisson")
    model.train(1, train_size=0.5)
    z1 = model.get_latent_representation(give_mean=True, mc_samples=1)
    assert z1.ndim == 2
    z2 = model.get_latent_representation(give_mean=False, mc_samples=15)
    assert (z2.ndim == 3) and (z2.shape[0] == 15)


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
    adata2 = adata.copy()
    model.get_elbo(adata2)
    assert adata.uns[_constants._SCVI_UUID_KEY] != adata2.uns[_constants._SCVI_UUID_KEY]

    # test mismatched categories raises ValueError
    adata2 = synthetic_iid()
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
                "{} already exists. Please provide an unexisting directory for saving.".format(
                    dir_path
                )
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

    def test_save_load_model(cls, adata, save_path, prefix=None, legacy=False):
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
        if legacy:
            legacy_save(
                model, save_path, overwrite=True, save_anndata=True, prefix=prefix
            )
        else:
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
        assert model.adata_manager.data_registry["batch"] == dict(
            attr_name="obs", attr_key="_scvi_batch", mod_key=None
        )

        z2 = model.get_latent_representation()
        test_idx2 = model.validation_indices
        np.testing.assert_array_equal(z1, z2)
        np.testing.assert_array_equal(test_idx1, test_idx2)
        assert model.is_trained is True

    save_path = os.path.join(save_path, "tmp")
    adata = synthetic_iid()

    for cls in [SCVI, LinearSCVI, TOTALVI, PEAKVI]:
        test_save_load_model(
            cls, adata, save_path, prefix=f"{cls.__name__}_", legacy=True
        )
        test_save_load_model(cls, adata, save_path, prefix=f"{cls.__name__}_")
        # Test load prioritizes newer save paradigm and thus mismatches legacy save.
        with pytest.raises(AssertionError):
            test_save_load_model(
                cls, adata, save_path, prefix=f"{cls.__name__}_", legacy=True
            )

    # AUTOZI
    def test_save_load_autozi(legacy=False):
        prefix = "AUTOZI_"
        model = AUTOZI(adata, latent_distribution="normal")
        model.train(1, train_size=0.5)
        ab1 = model.get_alphas_betas()
        if legacy:
            legacy_save(
                model, save_path, overwrite=True, save_anndata=True, prefix=prefix
            )
        else:
            model.save(save_path, overwrite=True, save_anndata=True, prefix=prefix)
            model.view_setup_args(save_path, prefix=prefix)
        model = AUTOZI.load(save_path, prefix=prefix)
        model.get_latent_representation()
        tmp_adata = scvi.data.synthetic_iid(n_genes=200)
        with pytest.raises(ValueError):
            AUTOZI.load(save_path, adata=tmp_adata, prefix=prefix)
        model = AUTOZI.load(save_path, adata=adata, prefix=prefix)
        assert "batch" in model.adata_manager.data_registry
        assert model.adata_manager.data_registry["batch"] == dict(
            attr_name="obs", attr_key="_scvi_batch", mod_key=None
        )

        ab2 = model.get_alphas_betas()
        np.testing.assert_array_equal(ab1["alpha_posterior"], ab2["alpha_posterior"])
        np.testing.assert_array_equal(ab1["beta_posterior"], ab2["beta_posterior"])
        assert model.is_trained is True

    AUTOZI.setup_anndata(adata, batch_key="batch", labels_key="labels")
    test_save_load_autozi(legacy=True)
    test_save_load_autozi()
    # Test load prioritizes newer save paradigm and thus mismatches legacy save.
    with pytest.raises(AssertionError):
        test_save_load_autozi(legacy=True)

    # SCANVI
    def test_save_load_scanvi(legacy=False):
        prefix = "SCANVI_"
        model = SCANVI(adata)
        model.train(max_epochs=1, train_size=0.5)
        p1 = model.predict()
        if legacy:
            legacy_save(
                model, save_path, overwrite=True, save_anndata=True, prefix=prefix
            )
        else:
            model.save(save_path, overwrite=True, save_anndata=True, prefix=prefix)
            model.view_setup_args(save_path, prefix=prefix)
        model = SCANVI.load(save_path, prefix=prefix)
        model.get_latent_representation()
        tmp_adata = scvi.data.synthetic_iid(n_genes=200)
        with pytest.raises(ValueError):
            SCANVI.load(save_path, adata=tmp_adata, prefix=prefix)
        model = SCANVI.load(save_path, adata=adata, prefix=prefix)
        assert "batch" in model.adata_manager.data_registry
        assert model.adata_manager.data_registry["batch"] == dict(
            attr_name="obs", attr_key="_scvi_batch", mod_key=None
        )

        p2 = model.predict()
        np.testing.assert_array_equal(p1, p2)
        assert model.is_trained is True

    SCANVI.setup_anndata(adata, "labels", "label_0", batch_key="batch")
    test_save_load_scanvi(legacy=True)
    test_save_load_scanvi()
    # Test load prioritizes newer save paradigm and thus mismatches legacy save.
    with pytest.raises(AssertionError):
        test_save_load_scanvi(legacy=True)


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
    adata3 = adata.copy()

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
    adata2_manager = manager_from_setup_dict(SCVI, adata2, LEGACY_SETUP_DICT)
    np.testing.assert_equal(
        field_registries_legacy_subset,
        adata2_manager.registry[_constants._FIELD_REGISTRIES_KEY],
    )

    # Test transfer.
    adata3_manager = adata_manager.transfer_fields(adata3)
    np.testing.assert_equal(
        field_registries,
        adata3_manager.registry[_constants._FIELD_REGISTRIES_KEY],
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
    a = scvi.data.synthetic_iid()
    m = scvi.model.SCVI.load(pretrained_scvi_path, adata=a)
    m.train(1)
    pretrained_totalvi_path = os.path.join(save_path, "testing_models/080_totalvi")
    m = scvi.model.TOTALVI.load(pretrained_totalvi_path, adata=a)
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


def test_ann_dataloader():
    a = scvi.data.synthetic_iid()
    adata_manager = generic_setup_adata_manager(
        a, batch_key="batch", labels_key="labels"
    )

    # test that batch sampler drops the last batch if it has less than 3 cells
    assert a.n_obs == 400
    adl = AnnDataLoader(adata_manager, batch_size=397, drop_last=3)
    assert len(adl) == 2
    for i, x in enumerate(adl):
        pass
    assert i == 1
    adl = AnnDataLoader(adata_manager, batch_size=398, drop_last=3)
    assert len(adl) == 1
    for i, x in enumerate(adl):
        pass
    assert i == 0
    with pytest.raises(ValueError):
        AnnDataLoader(adata_manager, batch_size=1, drop_last=2)


def test_semisupervised_dataloader():
    # test label resampling
    n_samples_per_label = 10
    a = synthetic_iid()
    adata_manager = generic_setup_adata_manager(
        a, batch_key="batch", labels_key="labels"
    )
    dl = SemiSupervisedDataLoader(
        adata_manager,
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
    ds = DeviceBackedDataSplitter(adata_manager, train_size=1.0, use_gpu=None)
    ds.setup()
    train_dl = ds.train_dataloader()
    ds.val_dataloader()
    loaded_x = next(iter(train_dl))["X"]
    assert len(loaded_x) == a.shape[0]
    np.testing.assert_array_equal(loaded_x.cpu().numpy(), a.X)

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
    adata_manager = generic_setup_adata_manager(
        a, batch_key="batch", labels_key="labels"
    )
    ds = SemiSupervisedDataSplitter(adata_manager, "asdf")
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
    ds = SemiSupervisedDataSplitter(adata_manager, unknown_label)
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
    SCANVI.setup_anndata(
        adata,
        "labels",
        "label_0",
        batch_key="batch",
    )
    model = SCANVI(adata, n_latent=10)
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
        autozivae.train(1, plan_kwargs=dict(lr=1e-2), check_val_every_n_epoch=1)
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
        autozivae.train(1, plan_kwargs=dict(lr=1e-2), check_val_every_n_epoch=1)
        assert hasattr(autozivae.module, "library_log_means") and hasattr(
            autozivae.module, "library_log_vars"
        )
        assert len(autozivae.history["elbo_train"]) == 1
        assert len(autozivae.history["elbo_validation"]) == 1
        autozivae.get_elbo(indices=autozivae.validation_indices)
        autozivae.get_reconstruction_error(indices=autozivae.validation_indices)
        autozivae.get_marginal_ll(indices=autozivae.validation_indices, n_mc_samples=3)
        autozivae.get_alphas_betas()


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
    TOTALVI.setup_anndata(
        adata,
        batch_key="batch",
        protein_expression_obsm_key="protein_expression",
        protein_names_uns_key="protein_names",
    )
    model = TOTALVI(adata)
    adata2 = synthetic_iid()
    model.get_elbo(adata2)

    # test that we catch incorrect mappings
    adata = synthetic_iid()
    TOTALVI.setup_anndata(
        adata,
        batch_key="batch",
        protein_expression_obsm_key="protein_expression",
        protein_names_uns_key="protein_names",
    )
    adata2 = synthetic_iid()
    adata2.obs.batch.cat.rename_categories(["batch_0", "batch_10"], inplace=True)
    with pytest.raises(ValueError):
        model.get_elbo(adata2)

    # test that same mapping different order is okay
    adata = synthetic_iid()
    TOTALVI.setup_anndata(
        adata,
        batch_key="batch",
        protein_expression_obsm_key="protein_expression",
        protein_names_uns_key="protein_names",
    )
    adata2 = synthetic_iid()
    adata2.obs.batch.cat.rename_categories(["batch_1", "batch_0"], inplace=True)
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
    CondSCVI.setup_anndata(dataset, labels_key="labels")
    sc_model = CondSCVI(dataset, n_latent=n_latent, n_layers=n_layers)
    sc_model.train(1, train_size=1)

    # step 2 learn destVI with multiple amortization scheme

    for amor_scheme in ["both", "none", "proportion", "latent"]:
        DestVI.setup_anndata(dataset, layer=None)
        spatial_model = DestVI.from_rna_model(
            dataset,
            sc_model,
            amortization=amor_scheme,
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


def test_early_stopping():
    n_epochs = 100

    adata = synthetic_iid()
    SCVI.setup_anndata(
        adata,
        batch_key="batch",
        labels_key="labels",
    )
    model = SCVI(adata)
    model.train(n_epochs, early_stopping=True, plan_kwargs=dict(lr=0))
    assert len(model.history["elbo_train"]) < n_epochs
