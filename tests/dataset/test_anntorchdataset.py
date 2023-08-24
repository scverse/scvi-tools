import numpy as np
import pytest
from torch.utils.data import Dataset

import scvi
from scvi import REGISTRY_KEYS
from scvi.utils._exceptions import InvalidParameterError
from tests.dataset.utils import generic_setup_adata_manager


def test_anntorchdataset_init():
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


def test_anntorchdataset_default_dtypes():
    adata = scvi.data.synthetic_iid()
    manager = generic_setup_adata_manager(adata, batch_key="batch")

    dataset = manager.create_torch_dataset()
    batch = dataset[:10]
    assert isinstance(batch, dict)
    assert REGISTRY_KEYS.X_KEY in batch
    assert REGISTRY_KEYS.BATCH_KEY in batch
    assert isinstance(batch[REGISTRY_KEYS.X_KEY], np.ndarray)
    assert isinstance(batch[REGISTRY_KEYS.BATCH_KEY], np.ndarray)
    assert batch[REGISTRY_KEYS.X_KEY].dtype == np.float32
    assert batch[REGISTRY_KEYS.BATCH_KEY].dtype == np.int64


def test_anntorchdataset_getitem_tensors():
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
    dataset = manager.create_torch_dataset(
        data_and_attributes={REGISTRY_KEYS.X_KEY: np.float64}
    )
    assert isinstance(dataset.keys_and_dtypes, dict)
    assert list(dataset.keys_and_dtypes.keys()) == [REGISTRY_KEYS.X_KEY]
    assert dataset.keys_and_dtypes[REGISTRY_KEYS.X_KEY] == np.float64

    with pytest.raises(KeyError):
        manager.create_torch_dataset(data_and_attributes=["not_a_key"])

    with pytest.raises(KeyError):
        manager.create_torch_dataset(data_and_attributes=[REGISTRY_KEYS.CAT_COVS_KEY])

    with pytest.raises(InvalidParameterError):
        manager.create_torch_dataset(data_and_attributes=1)
