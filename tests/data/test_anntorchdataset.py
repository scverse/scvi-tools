from __future__ import annotations

import os

import anndata
import numpy as np
import pytest
import torch
from torch.utils.data import Dataset

import scvi
from scvi import REGISTRY_KEYS
from tests.data.utils import generic_setup_adata_manager


def test_init():
    adata = scvi.data.synthetic_iid()
    manager = generic_setup_adata_manager(adata)

    dataset = manager.create_torch_dataset()
    assert dataset is not None
    assert isinstance(dataset, Dataset)
    assert len(dataset) > 0
    assert dataset.adata_manager is manager
    assert isinstance(dataset.keys_and_dtypes, dict)
    assert len(dataset.data) > 0
    assert len(dataset.registered_keys) > 0


def test_default_dtypes():
    adata = scvi.data.synthetic_iid()
    manager = generic_setup_adata_manager(adata, batch_key="batch")

    dataset = manager.create_torch_dataset()
    data = dataset[:10]
    assert isinstance(data, dict)
    assert REGISTRY_KEYS.X_KEY in data
    assert REGISTRY_KEYS.BATCH_KEY in data
    assert isinstance(data[REGISTRY_KEYS.X_KEY], np.ndarray)
    assert isinstance(data[REGISTRY_KEYS.BATCH_KEY], np.ndarray)
    assert data[REGISTRY_KEYS.X_KEY].dtype == np.float32
    assert data[REGISTRY_KEYS.BATCH_KEY].dtype == np.int64


def test_getitem_tensors():
    adata = scvi.data.synthetic_iid()
    manager = generic_setup_adata_manager(adata, batch_key="batch")

    # default
    dataset = manager.create_torch_dataset()
    assert isinstance(dataset.keys_and_dtypes, dict)
    assert REGISTRY_KEYS.X_KEY in dataset.keys_and_dtypes
    assert REGISTRY_KEYS.BATCH_KEY in dataset.keys_and_dtypes

    # list
    dataset = manager.create_torch_dataset(data_and_attributes=[REGISTRY_KEYS.X_KEY])
    assert isinstance(dataset.keys_and_dtypes, dict)
    assert list(dataset.keys_and_dtypes.keys()) == [REGISTRY_KEYS.X_KEY]

    # dict
    dataset = manager.create_torch_dataset(data_and_attributes={REGISTRY_KEYS.X_KEY: np.float64})
    assert isinstance(dataset.keys_and_dtypes, dict)
    assert list(dataset.keys_and_dtypes.keys()) == [REGISTRY_KEYS.X_KEY]
    assert dataset.keys_and_dtypes[REGISTRY_KEYS.X_KEY] == np.float64

    with pytest.raises(KeyError):
        manager.create_torch_dataset(data_and_attributes=["not_a_key"])

    with pytest.raises(KeyError):
        manager.create_torch_dataset(data_and_attributes=[REGISTRY_KEYS.CAT_COVS_KEY])

    with pytest.raises(ValueError):
        manager.create_torch_dataset(data_and_attributes=1)


def test_getitem(n_genes: int = 50):
    adata = scvi.data.synthetic_iid(n_genes=n_genes)
    manager = generic_setup_adata_manager(adata, batch_key="batch")

    dataset = manager.create_torch_dataset(load_sparse_tensor=True)

    data = dataset[:10]
    assert isinstance(data[REGISTRY_KEYS.X_KEY], np.ndarray)
    assert data[REGISTRY_KEYS.X_KEY].dtype == np.float32
    assert data[REGISTRY_KEYS.X_KEY].shape == (10, n_genes)
    assert isinstance(data[REGISTRY_KEYS.BATCH_KEY], np.ndarray)
    assert data[REGISTRY_KEYS.BATCH_KEY].dtype == np.int64
    assert data[REGISTRY_KEYS.BATCH_KEY].shape == (10, 1)

    data = dataset[[0, 1, 2]]
    assert isinstance(data[REGISTRY_KEYS.X_KEY], np.ndarray)
    assert data[REGISTRY_KEYS.X_KEY].dtype == np.float32
    assert data[REGISTRY_KEYS.X_KEY].shape == (3, n_genes)
    assert isinstance(data[REGISTRY_KEYS.BATCH_KEY], np.ndarray)
    assert data[REGISTRY_KEYS.BATCH_KEY].dtype == np.int64
    assert data[REGISTRY_KEYS.BATCH_KEY].shape == (3, 1)

    data = dataset[1]
    assert isinstance(data[REGISTRY_KEYS.X_KEY], np.ndarray)
    assert data[REGISTRY_KEYS.X_KEY].dtype == np.float32
    assert data[REGISTRY_KEYS.X_KEY].shape == (1, n_genes)
    assert isinstance(data[REGISTRY_KEYS.BATCH_KEY], np.ndarray)
    assert data[REGISTRY_KEYS.BATCH_KEY].dtype == np.int64
    assert data[REGISTRY_KEYS.BATCH_KEY].shape == (1, 1)


@pytest.mark.parametrize("sparse_format", ["csr_matrix", "csc_matrix"])
def test_load_sparse_tensor(sparse_format: str | None, n_genes: int = 50):
    adata = scvi.data.synthetic_iid(sparse_format=sparse_format, n_genes=n_genes)
    manager = generic_setup_adata_manager(adata, batch_key="batch")

    dataset = manager.create_torch_dataset(load_sparse_tensor=True)

    data = dataset[:10]
    assert isinstance(data[REGISTRY_KEYS.X_KEY], torch.Tensor)
    assert data[REGISTRY_KEYS.X_KEY].dtype == torch.float32
    assert data[REGISTRY_KEYS.X_KEY].shape == (10, n_genes)
    assert isinstance(data[REGISTRY_KEYS.BATCH_KEY], np.ndarray)
    assert data[REGISTRY_KEYS.BATCH_KEY].dtype == np.int64
    assert data[REGISTRY_KEYS.BATCH_KEY].shape == (10, 1)

    data = dataset[[0, 1, 2]]
    assert isinstance(data[REGISTRY_KEYS.X_KEY], torch.Tensor)
    assert data[REGISTRY_KEYS.X_KEY].dtype == torch.float32
    assert data[REGISTRY_KEYS.X_KEY].shape == (3, n_genes)
    assert isinstance(data[REGISTRY_KEYS.BATCH_KEY], np.ndarray)
    assert data[REGISTRY_KEYS.BATCH_KEY].dtype == np.int64
    assert data[REGISTRY_KEYS.BATCH_KEY].shape == (3, 1)

    data = dataset[1]
    assert isinstance(data[REGISTRY_KEYS.X_KEY], torch.Tensor)
    assert data[REGISTRY_KEYS.X_KEY].dtype == torch.float32
    assert data[REGISTRY_KEYS.X_KEY].shape == (1, n_genes)
    assert isinstance(data[REGISTRY_KEYS.BATCH_KEY], np.ndarray)
    assert data[REGISTRY_KEYS.BATCH_KEY].dtype == np.int64
    assert data[REGISTRY_KEYS.BATCH_KEY].shape == (1, 1)


def test_load_sparse_tensor_backed(save_path: str, n_genes: int = 50):
    adata = scvi.data.synthetic_iid(sparse_format="csr_matrix", n_genes=n_genes)
    adata_path = os.path.join(save_path, "adata.h5ad")
    adata.write(adata_path)
    del adata

    adata = anndata.read_h5ad(adata_path, backed="r")
    manager = generic_setup_adata_manager(adata, batch_key="batch")

    dataset = manager.create_torch_dataset(load_sparse_tensor=True)
    assert dataset.adata_manager.adata.isbacked

    data = dataset[:10]
    assert isinstance(data[REGISTRY_KEYS.X_KEY], torch.Tensor)
    assert data[REGISTRY_KEYS.X_KEY].dtype == torch.float32
    assert data[REGISTRY_KEYS.X_KEY].shape == (10, n_genes)
    assert isinstance(data[REGISTRY_KEYS.BATCH_KEY], np.ndarray)
    assert data[REGISTRY_KEYS.BATCH_KEY].dtype == np.int64
    assert data[REGISTRY_KEYS.BATCH_KEY].shape == (10, 1)

    data = dataset[[0, 1, 2]]
    assert isinstance(data[REGISTRY_KEYS.X_KEY], torch.Tensor)
    assert data[REGISTRY_KEYS.X_KEY].dtype == torch.float32
    assert data[REGISTRY_KEYS.X_KEY].shape == (3, n_genes)
    assert isinstance(data[REGISTRY_KEYS.BATCH_KEY], np.ndarray)
    assert data[REGISTRY_KEYS.BATCH_KEY].dtype == np.int64
    assert data[REGISTRY_KEYS.BATCH_KEY].shape == (3, 1)

    data = dataset[1]
    assert isinstance(data[REGISTRY_KEYS.X_KEY], torch.Tensor)
    assert data[REGISTRY_KEYS.X_KEY].dtype == torch.float32
    assert data[REGISTRY_KEYS.X_KEY].shape == (1, n_genes)
    assert isinstance(data[REGISTRY_KEYS.BATCH_KEY], np.ndarray)
    assert data[REGISTRY_KEYS.BATCH_KEY].dtype == np.int64
    assert data[REGISTRY_KEYS.BATCH_KEY].shape == (1, 1)
