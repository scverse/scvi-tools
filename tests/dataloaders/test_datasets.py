import os

import anndata
import mudata
import numpy as np
import pytest
import torch

from scvi import REGISTRY_KEYS
from scvi.data import synthetic_iid
from scvi.dataloaders import AnnTorchDataset
from tests.dataset.utils import (
    generic_setup_adata_manager,
    generic_setup_mudata_manager,
)


def test_init_anntorchdataset():
    adata = synthetic_iid()
    manager = generic_setup_adata_manager(adata, batch_key="batch", labels_key="labels")

    dataset = AnnTorchDataset(manager)
    assert len(dataset) == adata.n_obs

    ########################################
    assert hasattr(dataset, "data")
    assert isinstance(dataset.data, dict)
    assert all([isinstance(key, str) for key in dataset.data.keys()])
    assert all([isinstance(val, np.ndarray) for val in dataset.data.values()])
    ########################################

    ########################################
    assert hasattr(dataset, "attributes_and_types")
    assert isinstance(dataset.attributes_and_types, dict)
    assert all([isinstance(key, str) for key in dataset.attributes_and_types.keys()])
    assert all([isinstance(val, type) for val in dataset.attributes_and_types.values()])
    ########################################


def test_getitem_tensors_anntorchdataset():
    adata = synthetic_iid()
    manager = generic_setup_adata_manager(adata, batch_key="batch", labels_key="labels")

    ########################################
    # default getitem_tensors
    getitem_tensors = None
    dataset = AnnTorchDataset(manager, getitem_tensors=getitem_tensors)
    expected_keys = [
        REGISTRY_KEYS.X_KEY,
        REGISTRY_KEYS.BATCH_KEY,
        REGISTRY_KEYS.LABELS_KEY,
    ]
    data = dataset[0]
    assert isinstance(data, dict)
    assert all([(key in expected_keys) for key in data.keys()])
    ########################################

    ########################################
    # getitem_tensors list
    getitem_tensors = [REGISTRY_KEYS.X_KEY, REGISTRY_KEYS.BATCH_KEY]
    dataset = AnnTorchDataset(manager, getitem_tensors=getitem_tensors)
    expected_keys = [REGISTRY_KEYS.X_KEY, REGISTRY_KEYS.BATCH_KEY]
    data = dataset[0]
    assert isinstance(data, dict)
    assert all([(key in expected_keys) for key in data.keys()])
    X, batch = data[REGISTRY_KEYS.X_KEY], data[REGISTRY_KEYS.BATCH_KEY]
    assert isinstance(X, np.ndarray)
    assert isinstance(batch, np.ndarray)
    assert X.dtype == np.float32
    assert batch.dtype == np.float32
    ## invalid key
    getitem_tensors = [REGISTRY_KEYS.X_KEY, "invalid_key"]
    with pytest.raises(ValueError):
        _ = AnnTorchDataset(manager, getitem_tensors=getitem_tensors)
    ########################################

    ########################################
    # getitem_tensors dict
    getitem_tensors = {
        REGISTRY_KEYS.X_KEY: np.float32,
        REGISTRY_KEYS.BATCH_KEY: np.int64,
    }
    dataset = AnnTorchDataset(manager, getitem_tensors=getitem_tensors)
    expected_keys = [REGISTRY_KEYS.X_KEY, REGISTRY_KEYS.BATCH_KEY]
    data = dataset[0]
    assert isinstance(data, dict)
    assert all([(key in expected_keys) for key in data.keys()])
    X, batch = data[REGISTRY_KEYS.X_KEY], data[REGISTRY_KEYS.BATCH_KEY]
    assert isinstance(X, np.ndarray)
    assert isinstance(batch, np.ndarray)
    assert X.dtype == np.float32
    assert batch.dtype == np.int64
    ## invalid key
    getitem_tensors = {
        REGISTRY_KEYS.X_KEY: np.float32,
        "invalid_key": np.int64,
    }
    with pytest.raises(ValueError):
        _ = AnnTorchDataset(manager, getitem_tensors=getitem_tensors)
    ## invalid type
    getitem_tensors = {
        REGISTRY_KEYS.X_KEY: np.float32,
        REGISTRY_KEYS.BATCH_KEY: "invalid_type",
    }
    with pytest.raises(ValueError):
        _ = AnnTorchDataset(manager, getitem_tensors=getitem_tensors)
    ########################################


@pytest.mark.parametrize("sparse", [True, False])
def test_getitem_anntorchdataset(
    sparse: bool, batch_size: int = 64, n_batches: int = 2, n_genes: int = 25
):
    adata = synthetic_iid(
        batch_size=batch_size, n_batches=n_batches, n_genes=n_genes, sparse=sparse
    )
    manager = generic_setup_adata_manager(adata, batch_key="batch")

    dataset = AnnTorchDataset(manager)

    def _check_data(data: dict):
        X, batch = data[REGISTRY_KEYS.X_KEY], data[REGISTRY_KEYS.BATCH_KEY]
        assert X.shape[0] == batch.shape[0]
        assert X.ndim == 2
        assert batch.ndim == 2
        assert X.shape[1] == n_genes
        assert batch.shape[1] == 1

    ########################################
    # getitem single integer index
    for i in range(len(dataset)):
        _check_data(dataset[i])
    for i in range(-len(dataset), 0):
        _check_data(dataset[i])
    ## invalid index
    with pytest.raises(IndexError):
        _ = dataset[len(dataset)]
    ########################################

    ########################################
    # getitem list
    indices = [0, 1, 2]
    _check_data(dataset[indices])
    indices = [3, 2, 1]
    _check_data(dataset[indices])
    ## invalid index
    indices = [0, 1, 2, len(dataset)]
    with pytest.raises(IndexError):
        _ = dataset[indices]
    ########################################

    ########################################
    # getitem np.ndarray
    indices = np.array([0, 1, 2])
    _check_data(dataset[indices])
    indices = np.array([3, 2, 1])
    _check_data(dataset[indices])
    ## invalid index
    indices = np.array([0, 1, 2, len(dataset)])
    with pytest.raises(IndexError):
        _ = dataset[indices]

    ########################################
    # getitem slice
    _check_data(dataset[:3])
    _check_data(dataset[3:])
    _check_data(dataset[3:10])
    _check_data(dataset[3:10:2])
    ########################################


@pytest.mark.parametrize("sparse", [True, False])
def test_disk_backed_anntorchdataset(save_path: str, sparse: bool):
    adata = synthetic_iid(sparse=sparse)
    adata_path = os.path.join(
        save_path, f"disk_backed_anntorchdataset_adata_{sparse}.h5ad"
    )
    adata.write(adata_path)
    del adata

    adata = anndata.read_h5ad(adata_path, backed="r")
    manager = generic_setup_adata_manager(adata, batch_key="batch", labels_key="labels")

    dataset = AnnTorchDataset(manager)
    expected_keys = [
        REGISTRY_KEYS.X_KEY,
        REGISTRY_KEYS.BATCH_KEY,
        REGISTRY_KEYS.LABELS_KEY,
    ]

    ########################################
    assert adata.isbacked
    assert hasattr(dataset, "data")
    assert isinstance(dataset.data, dict)
    assert all([key in expected_keys for key in dataset.data.keys()])
    ########################################


@pytest.mark.parametrize("sparse", [True, False])
def test_cuda_backed_anntorchdataset(cuda: bool, sparse: bool):
    adata = synthetic_iid(sparse=sparse)
    manager = generic_setup_adata_manager(adata, batch_key="batch", labels_key="labels")

    dataset = AnnTorchDataset(manager, accelerator="cuda", device_backed=True)
    expected_keys = [
        REGISTRY_KEYS.X_KEY,
        REGISTRY_KEYS.BATCH_KEY,
        REGISTRY_KEYS.LABELS_KEY,
    ]

    ########################################
    assert hasattr(dataset, "data")
    assert dataset.device.type == "cuda"
    assert isinstance(dataset.data, dict)
    assert all([key in expected_keys for key in dataset.data.keys()])
    assert all([isinstance(value, torch.Tensor) for value in dataset.data.values()])
    assert all([value.device.type == "cuda" for value in dataset.data.values()])
    ########################################


def test_obsm_anntorchdataset(n_proteins: int = 20):
    adata = synthetic_iid(n_proteins=n_proteins)
    manager = generic_setup_adata_manager(
        adata,
        batch_key="batch",
        labels_key="labels",
        protein_expression_obsm_key="protein_expression",
    )

    dataset = AnnTorchDataset(manager)
    expected_keys = [
        REGISTRY_KEYS.X_KEY,
        REGISTRY_KEYS.BATCH_KEY,
        REGISTRY_KEYS.LABELS_KEY,
        REGISTRY_KEYS.PROTEIN_EXP_KEY,
    ]

    ########################################
    assert all([key in expected_keys for key in dataset.data.keys()])
    assert dataset.data[REGISTRY_KEYS.PROTEIN_EXP_KEY].shape[0] == adata.n_obs
    assert dataset.data[REGISTRY_KEYS.PROTEIN_EXP_KEY].shape[1] == n_proteins
    ########################################


def test_mudata_anntorchdataset():
    adata = synthetic_iid()
    bdata = synthetic_iid()
    mdata = mudata.MuData({"rna": adata, "protein": bdata})
    manager = generic_setup_mudata_manager(
        mdata,
        layer_mod="rna",
        layer=None,
        batch_mod="rna",
        batch_key="batch",
        protein_expression_mod="protein",
        protein_expression_layer=None,
    )

    dataset = AnnTorchDataset(manager)
    expected_keys = [
        REGISTRY_KEYS.X_KEY,
        REGISTRY_KEYS.BATCH_KEY,
        REGISTRY_KEYS.PROTEIN_EXP_KEY,
    ]

    ########################################
    assert all([key in expected_keys for key in dataset.data.keys()])
    assert dataset.data[REGISTRY_KEYS.PROTEIN_EXP_KEY].shape[0] == adata.n_obs
    assert dataset.data[REGISTRY_KEYS.PROTEIN_EXP_KEY].shape[1] == bdata.n_vars
    ########################################
