import numpy as np
import pytest

import scvi
from scvi.data import synthetic_iid
from scvi.data._compat import LEGACY_REGISTRY_KEY_MAP
from scvi.dataloaders import (
    AnnDataLoader,
    DataSplitter,
    SemiSupervisedDataLoader,
    SemiSupervisedDataSplitter,
)
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


def test_ann_dataloader():
    a = scvi.data.synthetic_iid()
    adata_manager = generic_setup_adata_manager(
        a, batch_key="batch", labels_key="labels"
    )

    # test that batch sampler drops the last batch if it has less than 3 cells
    assert a.n_obs == 400
    adl = AnnDataLoader(adata_manager, batch_size=397, drop_last=3)
    assert len(adl) == 2
    for _i, _ in enumerate(adl):
        pass
    assert _i == 1
    adl = AnnDataLoader(adata_manager, batch_size=398, drop_last=3)
    assert len(adl) == 1
    for _i, _ in enumerate(adl):
        pass
    assert _i == 0
    with pytest.raises(ValueError):
        AnnDataLoader(adata_manager, batch_size=1, drop_last=2)


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


def test_default_data_loader_kwargs():
    a = synthetic_iid()
    adata_manager = generic_setup_adata_manager(
        a, batch_key="batch", labels_key="labels"
    )
    # test that data_loader_kwargs get default drop_last=3 when respective keyword arg is not provided to DataSplitter
    ds = DataSplitter(adata_manager, train_size=400 / a.n_obs, batch_size=397)
    ds.setup()
    adl = ds.train_dataloader()
    assert len(adl) == 2
    for i, _x in enumerate(adl):
        pass
    assert i == 1

    ds = DataSplitter(adata_manager, train_size=400 / a.n_obs, batch_size=398)
    ds.setup()
    adl = ds.train_dataloader()
    assert len(adl) == 1
    for i, _x in enumerate(adl):
        pass
    assert i == 0

    # test that splitter kwargs are able to override default values for data_loader_kwargs: integer drop last batch size
    ds = DataSplitter(
        adata_manager, train_size=400 / a.n_obs, batch_size=398, drop_last=2
    )
    ds.setup()
    adl = ds.train_dataloader()
    assert len(adl) == 2
    for i, _x in enumerate(adl):
        pass
    assert i == 1

    # test that splitter kwargs are able to override default values for data_loader_kwargs: boolean drop_last=True
    ds = DataSplitter(adata_manager, train_size=400 / a.n_obs, batch_size=398)
    ds.setup()
    adl = ds.train_dataloader()
    assert len(adl) == 1
    for i, _x in enumerate(adl):
        pass
    assert i == 0

    with pytest.raises(ValueError):
        ds = DataSplitter(
            adata_manager, train_size=400 / a.n_obs, batch_size=1, drop_last=2
        )
        ds.setup()
        adl = ds.train_dataloader()


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
