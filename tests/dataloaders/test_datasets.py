from itertools import product

import numpy as np
import pytest

from scvi.data import synthetic_iid
from scvi.dataloaders import AnnTorchDataset
from tests.dataset.utils import generic_setup_adata_manager


@pytest.mark.parametrize(
    "batch_size, n_batches, sparse", product([128, 256], [1, 2], [True, False])
)
def test_basic_anntorchdataset(batch_size: int, n_batches: int, sparse: bool):
    adata = synthetic_iid(batch_size=batch_size, n_batches=n_batches, sparse=sparse)
    manager = generic_setup_adata_manager(adata, batch_key="batch", labels_key="labels")

    # test init
    dataset = AnnTorchDataset(manager)
    assert hasattr(dataset, "data")  # data should be loaded on init
    assert len(dataset) == n_batches * batch_size
    attr_and_types = {k: np.float32 for k in manager.data_registry}
    assert dataset.attributes_and_types == attr_and_types

    # test getitem_tensors with list
    dataset = AnnTorchDataset(manager, getitem_tensors=["X", "batch"])
    assert dataset.attributes_and_types == {"X": np.float32, "batch": np.float32}
    ## test getitem with a single index
    data = dataset[0]
    X, batch = data["X"], data["batch"]
    assert X.dtype == np.float32
    assert batch.dtype == np.float32
    assert X.shape == (1, 100)  # single observation should be a row vector
    assert batch.shape == (1, 1)
    ## test getitem with a slice
    data = dataset[:10]
    X, batch = data["X"], data["batch"]
    assert X.shape == (10, 100)
    assert batch.shape == (10, 1)

    # test getitem_tensors with dict
    getitem_tensors = {"X": np.float64, "batch": np.int32}
    dataset = AnnTorchDataset(manager, getitem_tensors=getitem_tensors)
    assert dataset.attributes_and_types == getitem_tensors
    ## test getitem with a single index
    data = dataset[0]
    X, batch = data["X"], data["batch"]
    assert X.dtype == np.float64
    assert batch.dtype == np.int32
    assert X.shape == (1, 100)  # single observation should be a row vector
    assert batch.shape == (1, 1)
    ## test getitem with a slice
    data = dataset[:10]
    X, batch = data["X"], data["batch"]
    assert X.shape == (10, 100)
    assert batch.shape == (10, 1)


@pytest.mark.parametrize(
    "batch_size, n_batches, sparse", product([128, 256], [1, 2], [True, False])
)
def test_backed_adata_anntorchdataset(batch_size: int, n_batches: int, sparse: bool):
    pass


@pytest.mark.cuda
@pytest.mark.parametrize(
    "batch_size, n_batches, sparse", product([128, 256], [1, 2], [True, False])
)
def test_cuda_anntorchdataset(
    cuda: bool, batch_size: int, n_batches: int, sparse: bool
):
    pass
