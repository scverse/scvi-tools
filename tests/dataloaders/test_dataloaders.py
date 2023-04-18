from itertools import product

import numpy as np
import pytest
import torch
from torch.utils.data import SequentialSampler

from scvi import REGISTRY_KEYS
from scvi.data import synthetic_iid
from scvi.dataloaders import AnnDataLoader, AnnTorchDataset
from tests.dataset.utils import generic_setup_adata_manager


@pytest.mark.parametrize(
    "batch_size, n_batches, n_genes, sparse",
    product([128], [1, 2], [50], [True, False]),
)
def test_basic_anndataloader(
    batch_size: int, n_batches: int, n_genes: int, sparse: bool
):
    adata = synthetic_iid(
        batch_size=batch_size, n_batches=n_batches, n_genes=n_genes, sparse=sparse
    )
    manager = generic_setup_adata_manager(adata, batch_key="batch", labels_key="labels")

    dl = AnnDataLoader(manager, batch_size=batch_size)
    assert len(dl) == n_batches
    batch = next(iter(dl))
    for key, value in batch.items():
        assert key in (
            REGISTRY_KEYS.X_KEY,
            REGISTRY_KEYS.BATCH_KEY,
            REGISTRY_KEYS.LABELS_KEY,
        )
        assert isinstance(value, torch.Tensor)
        assert value.shape[0] == batch_size
        assert value.shape[1] == n_genes if key == REGISTRY_KEYS.X_KEY else 1


@pytest.mark.parametrize(
    "batch_size, n_batches, n_genes, sparse",
    product([128], [1, 2], [50], [True, False]),
)
def test_drop_last_anndataloader(
    batch_size: int, n_batches: int, n_genes: int, sparse: bool
):
    adata = synthetic_iid(
        batch_size=batch_size, n_batches=n_batches, n_genes=n_genes, sparse=sparse
    )
    manager = generic_setup_adata_manager(adata, batch_key="batch", labels_key="labels")

    # test drop_last with default sampler
    dl = AnnDataLoader(manager, batch_size=batch_size + 1, drop_last=True)
    assert len(dl) == (n_batches - 1)
    dl = AnnDataLoader(manager, batch_size=batch_size + 1, drop_last=False)
    assert len(dl) == n_batches
    for i, batch in enumerate(iter(dl)):
        if i < n_batches - 1:
            assert batch[REGISTRY_KEYS.X_KEY].shape[0] == batch_size + 1
        else:
            assert batch[REGISTRY_KEYS.X_KEY].shape[0] == batch_size - (n_batches - 1)

    # test drop_last with custom sampler
    sampler = SequentialSampler(AnnTorchDataset(manager))
    dl = AnnDataLoader(
        manager,
        sampler=sampler,
        batch_size=batch_size + 1,
        drop_last=True,
    )
    assert len(dl) == (n_batches - 1)
    dl = AnnDataLoader(
        manager,
        sampler=sampler,
        batch_size=batch_size + 1,
        drop_last=False,
    )
    assert len(dl) == n_batches


@pytest.mark.parametrize(
    "batch_size, n_batches, n_genes, sparse",
    product([128], [2], [50], [True, False]),
)
def test_iter_ndarray_anndataloader(
    batch_size: int, n_batches: int, n_genes: int, sparse: bool
):
    adata = synthetic_iid(
        batch_size=batch_size, n_batches=n_batches, n_genes=n_genes, sparse=sparse
    )
    manager = generic_setup_adata_manager(adata, batch_key="batch", labels_key="labels")

    dl = AnnDataLoader(manager, batch_size=batch_size, iter_ndarray=True)
    assert len(dl) == n_batches
    batch = next(iter(dl))
    for value in batch.values():
        assert isinstance(value, np.ndarray)

    with pytest.raises(ValueError):
        _ = AnnDataLoader(manager, iter_ndarray=True, device_backed=True)

    with pytest.raises(ValueError):
        _ = AnnDataLoader(manager, iter_ndarray=True, collate_fn=lambda x: x)


@pytest.mark.cuda
@pytest.mark.parametrize(
    "batch_size, n_batches, n_genes, sparse",
    product([128], [1, 2], [50], [True, False]),
)
def test_cuda_backed_anndataloader(
    cuda: bool, batch_size: int, n_batches: int, n_genes: int, sparse: bool
):
    assert cuda

    adata = synthetic_iid(
        batch_size=batch_size, n_batches=n_batches, n_genes=n_genes, sparse=sparse
    )
    manager = generic_setup_adata_manager(adata, batch_key="batch", labels_key="labels")

    dl = AnnDataLoader(
        manager, batch_size=batch_size, accelerator="cuda", device_backed=True
    )
    assert len(dl) == n_batches
    batch = next(iter(dl))
    for key, value in batch.items():
        assert key in (
            REGISTRY_KEYS.X_KEY,
            REGISTRY_KEYS.BATCH_KEY,
            REGISTRY_KEYS.LABELS_KEY,
        )
        assert value.shape[0] == batch_size
        assert value.shape[1] == n_genes if key == REGISTRY_KEYS.X_KEY else 1
        assert isinstance(value, torch.Tensor)
        assert value.device.type == "cuda"
