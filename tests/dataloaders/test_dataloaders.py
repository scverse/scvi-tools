import os
from math import ceil, floor

import anndata
import numpy as np
import pytest
import torch
from torch.utils.data import BatchSampler, RandomSampler, SequentialSampler

from scvi import REGISTRY_KEYS
from scvi.data import synthetic_iid
from scvi.dataloaders import (
    AnnDataLoader,
    AnnTorchDataset,
    ConcatDataLoader,
)
from tests.dataset.utils import generic_setup_adata_manager


def test_init_anndataloader():
    adata = synthetic_iid()
    manager = generic_setup_adata_manager(adata, batch_key="batch", labels_key="labels")

    ########################################
    dataloader = AnnDataLoader(manager)
    assert len(dataloader) == ceil(adata.n_obs / 128)
    ########################################

    ########################################
    dataset = AnnTorchDataset(manager)
    dataloader = AnnDataLoader(dataset=dataset)
    assert len(dataloader) == ceil(adata.n_obs / 128)
    ########################################

    ########################################
    with pytest.raises(ValueError):
        _ = AnnDataLoader(manager, dataset=dataset)

    with pytest.raises(ValueError):
        _ = AnnDataLoader()
    ########################################


def test_sampler_anndataloader():
    adata = synthetic_iid()
    manager = generic_setup_adata_manager(adata, batch_key="batch", labels_key="labels")
    dataset = AnnTorchDataset(manager)

    ########################################
    # default sampler
    dataloader = AnnDataLoader(manager, shuffle=True, batch_size=128, drop_last=True)
    assert dataloader.batch_size is None
    assert not dataloader.drop_last
    assert dataloader.batch_sampler is None
    assert isinstance(dataloader.sampler, BatchSampler)
    assert isinstance(dataloader.sampler.sampler, RandomSampler)

    dataloader = AnnDataLoader(manager, shuffle=False, batch_size=128, drop_last=False)
    assert dataloader.batch_size is None
    assert not dataloader.drop_last
    assert dataloader.batch_sampler is None
    assert isinstance(dataloader.sampler, BatchSampler)
    assert isinstance(dataloader.sampler.sampler, SequentialSampler)
    ########################################

    ########################################
    # custom sampler
    dataset = AnnTorchDataset(manager)
    sampler = SequentialSampler(dataset)
    dataloader = AnnDataLoader(
        dataset=dataset, sampler=sampler, shuffle=False, batch_size=10, drop_last=True
    )
    assert dataloader.batch_size == 10
    assert dataloader.drop_last
    assert isinstance(dataloader.sampler, SequentialSampler)
    assert isinstance(dataloader.batch_sampler, BatchSampler)
    assert dataloader.batch_sampler.batch_size == 10
    ########################################

    ########################################
    # sampler and shuffle are mutually exclusive
    sampler = RandomSampler(dataset)
    with pytest.raises(ValueError):
        _ = AnnDataLoader(
            dataset=dataset,
            sampler=sampler,
            shuffle=True,
            batch_size=10,
            drop_last=True,
        )
    ########################################

    ########################################
    # custom batch sampler
    sampler = RandomSampler(dataset)
    sampler = BatchSampler(sampler=sampler, batch_size=10, drop_last=True)
    dataloader = AnnDataLoader(dataset=dataset, sampler=sampler)
    assert dataloader.batch_size is None
    assert not dataloader.drop_last
    assert dataloader.batch_sampler is None
    assert isinstance(dataloader.sampler, BatchSampler)
    assert isinstance(dataloader.sampler.sampler, RandomSampler)
    assert dataloader.sampler.batch_size == 10
    ########################################


def test_indices_anndataloader():
    adata = synthetic_iid()
    manager = generic_setup_adata_manager(adata, batch_key="batch", labels_key="labels")

    ########################################
    indices = [0, 1, 2]
    dataloader = AnnDataLoader(manager, indices=indices)
    assert len(dataloader.dataset) == len(indices)
    ########################################

    ########################################
    indices = np.array([0, 1, 2])
    dataloader = AnnDataLoader(manager, indices=indices)
    assert len(dataloader.dataset) == len(indices)
    ########################################

    ########################################
    indices = np.arange(adata.n_obs // 2)
    dataloader = AnnDataLoader(manager, indices=indices)
    assert len(dataloader.dataset) == adata.n_obs // 2
    ########################################

    ########################################
    indices = np.array([False] * (adata.n_obs // 2) + [True] * (adata.n_obs // 2))
    dataloader = AnnDataLoader(manager, indices=indices)
    assert len(dataloader.dataset) == adata.n_obs // 2
    ########################################

    ########################################
    indices = np.arange(adata.n_obs + 1)
    with pytest.raises(ValueError):
        _ = AnnDataLoader(manager, indices=indices)


@pytest.mark.parametrize("sparse", [True, False])
def test_iter_anndataloader(
    sparse: bool, batch_size: int = 64, n_batches: int = 2, n_genes: int = 25
):
    adata = synthetic_iid(
        batch_size=batch_size, n_batches=n_batches, n_genes=n_genes, sparse=sparse
    )
    manager = generic_setup_adata_manager(adata, batch_key="batch", labels_key="labels")

    ########################################
    dataloader = AnnDataLoader(manager, batch_size=batch_size, drop_last=True)
    for batch in dataloader:
        X, batch, labels = (
            batch[REGISTRY_KEYS.X_KEY],
            batch[REGISTRY_KEYS.BATCH_KEY],
            batch[REGISTRY_KEYS.LABELS_KEY],
        )
        assert isinstance(X, torch.Tensor)
        assert isinstance(batch, torch.Tensor)
        assert isinstance(labels, torch.Tensor)
        assert X.shape == (batch_size, n_genes)
        assert batch.shape == (batch_size, 1)
        assert labels.shape == (batch_size, 1)
    ########################################

    ########################################
    dataloader = AnnDataLoader(manager, batch_size=adata.n_obs + 1, drop_last=True)
    with pytest.raises(StopIteration):
        _ = next(iter(dataloader))
    ########################################


@pytest.mark.parametrize("sparse", [True, False])
def test_iter_ndarray_anndataloader(
    sparse: bool, batch_size: int = 64, n_batches: int = 2, n_genes: int = 25
):
    adata = synthetic_iid(
        batch_size=batch_size, n_batches=n_batches, n_genes=n_genes, sparse=sparse
    )
    manager = generic_setup_adata_manager(adata, batch_key="batch", labels_key="labels")

    ########################################
    dataloader = AnnDataLoader(
        manager, batch_size=batch_size, drop_last=True, iter_ndarray=True
    )
    for batch in dataloader:
        X, batch, labels = (
            batch[REGISTRY_KEYS.X_KEY],
            batch[REGISTRY_KEYS.BATCH_KEY],
            batch[REGISTRY_KEYS.LABELS_KEY],
        )
        assert isinstance(X, np.ndarray)
        assert isinstance(batch, np.ndarray)
        assert isinstance(labels, np.ndarray)
        assert X.shape == (batch_size, n_genes)
        assert batch.shape == (batch_size, 1)
        assert labels.shape == (batch_size, 1)
    ########################################

    ########################################
    dataloader = AnnDataLoader(manager, batch_size=adata.n_obs + 1, drop_last=True)
    with pytest.raises(StopIteration):
        _ = next(iter(dataloader))
    ########################################


@pytest.mark.parametrize("sparse", [True, False])
def test_cuda_backed_anndataloader(
    cuda: bool,
    sparse: bool,
    batch_size: int = 64,
    n_batches: int = 2,
    n_genes: int = 25,
):
    assert cuda
    adata = synthetic_iid(
        batch_size=batch_size, n_batches=n_batches, n_genes=n_genes, sparse=sparse
    )
    manager = generic_setup_adata_manager(adata, batch_key="batch", labels_key="labels")

    ########################################
    dataloader = AnnDataLoader(
        manager,
        batch_size=batch_size,
        accelerator="cuda",
        device_backed=True,
        drop_last=True,
    )
    for batch in dataloader:
        X, batch, labels = (
            batch[REGISTRY_KEYS.X_KEY],
            batch[REGISTRY_KEYS.BATCH_KEY],
            batch[REGISTRY_KEYS.LABELS_KEY],
        )
        assert isinstance(X, torch.Tensor)
        assert isinstance(batch, torch.Tensor)
        assert isinstance(labels, torch.Tensor)
        assert X.device.type == "cuda"
        assert batch.device.type == "cuda"
        assert labels.device.type == "cuda"
        assert X.shape == (batch_size, n_genes)
        assert batch.shape == (batch_size, 1)
        assert labels.shape == (batch_size, 1)
    ########################################

    ########################################
    dataloader = AnnDataLoader(
        manager, batch_size=batch_size, drop_last=True, iter_ndarray=True
    )
    for batch in dataloader:
        X, batch, labels = (
            batch[REGISTRY_KEYS.X_KEY],
            batch[REGISTRY_KEYS.BATCH_KEY],
            batch[REGISTRY_KEYS.LABELS_KEY],
        )
        assert isinstance(X, np.ndarray)
        assert isinstance(batch, np.ndarray)
        assert isinstance(labels, np.ndarray)
        assert X.shape == (batch_size, n_genes)
        assert batch.shape == (batch_size, 1)
        assert labels.shape == (batch_size, 1)
    ########################################


@pytest.mark.parametrize("sparse", [True, False])
def test_disk_backed_anndataloader(
    save_path: str,
    sparse: bool,
    batch_size: int = 64,
    n_batches: int = 2,
    n_genes: int = 25,
):
    adata = synthetic_iid(
        batch_size=batch_size, n_batches=n_batches, n_genes=n_genes, sparse=sparse
    )
    adata_path = os.path.join(
        save_path, f"disk_backed_anndataloader_adata_{sparse}.h5ad"
    )
    adata.write(adata_path)
    del adata

    adata = anndata.read_h5ad(adata_path, backed="r")
    manager = generic_setup_adata_manager(adata, batch_key="batch", labels_key="labels")

    ########################################
    dataloader = AnnDataLoader(
        manager, batch_size=batch_size, accelerator="cuda", drop_last=True
    )
    for batch in dataloader:
        X, batch, labels = (
            batch[REGISTRY_KEYS.X_KEY],
            batch[REGISTRY_KEYS.BATCH_KEY],
            batch[REGISTRY_KEYS.LABELS_KEY],
        )
        assert isinstance(X, torch.Tensor)
        assert isinstance(batch, torch.Tensor)
        assert isinstance(labels, torch.Tensor)
        assert X.shape == (batch_size, n_genes)
        assert batch.shape == (batch_size, 1)
        assert labels.shape == (batch_size, 1)
    ########################################

    ########################################
    dataloader = AnnDataLoader(
        manager, batch_size=batch_size, drop_last=True, iter_ndarray=True
    )
    for batch in dataloader:
        X, batch, labels = (
            batch[REGISTRY_KEYS.X_KEY],
            batch[REGISTRY_KEYS.BATCH_KEY],
            batch[REGISTRY_KEYS.LABELS_KEY],
        )
        assert isinstance(X, np.ndarray)
        assert isinstance(batch, np.ndarray)
        assert isinstance(labels, np.ndarray)
        assert X.shape == (batch_size, n_genes)
        assert batch.shape == (batch_size, 1)
        assert labels.shape == (batch_size, 1)
    ########################################


@pytest.mark.parametrize("sparse", [True, False])
def test_init_concatdataloader(
    sparse: bool, batch_size: int = 64, n_batches: int = 2, n_genes: int = 25
):
    adata = synthetic_iid(
        batch_size=batch_size, n_batches=n_batches, n_genes=n_genes, sparse=sparse
    )
    manager = generic_setup_adata_manager(adata, batch_key="batch", labels_key="labels")

    ########################################
    indices_list = [
        np.arange(adata.n_obs // 2),
        np.arange(adata.n_obs // 2, adata.n_obs),
    ]
    dataloader = ConcatDataLoader(manager, indices_list=indices_list)
    assert len(dataloader) == ceil(adata.n_obs / 128)
    assert hasattr(dataloader, "indices_list")
    assert isinstance(dataloader.indices_list, list)
    assert len(dataloader.indices_list) == 2
    assert hasattr(dataloader, "dataloaders")
    assert isinstance(dataloader.dataloaders, list)
    assert all([isinstance(dl, AnnDataLoader) for dl in dataloader.dataloaders])
    ########################################


@pytest.mark.parametrize("sparse", [True, False])
def test_iter_concatdataloader(
    sparse: bool, batch_size: int = 64, n_batches: int = 2, n_genes: int = 25
):
    adata = synthetic_iid(
        batch_size=batch_size, n_batches=n_batches, n_genes=n_genes, sparse=sparse
    )
    manager = generic_setup_adata_manager(adata, batch_key="batch", labels_key="labels")

    ########################################
    # one list of indices
    indices_list = [np.arange(0, adata.n_obs // 2)]
    dataloader = ConcatDataLoader(
        manager, indices_list=indices_list, batch_size=10, drop_last=False
    )
    assert len(dataloader) == ceil((adata.n_obs / 2) / 10)
    dataloader = ConcatDataLoader(
        manager, indices_list=indices_list, batch_size=10, drop_last=True
    )
    assert len(dataloader) == floor((adata.n_obs / 2) / 10)
    for batch in dataloader:
        assert len(batch) == len(indices_list)
        X, batch, labels = (
            batch[0][REGISTRY_KEYS.X_KEY],
            batch[0][REGISTRY_KEYS.BATCH_KEY],
            batch[0][REGISTRY_KEYS.LABELS_KEY],
        )
        assert isinstance(X, torch.Tensor)
        assert isinstance(batch, torch.Tensor)
        assert isinstance(labels, torch.Tensor)
        assert X.shape == (10, n_genes)
        assert batch.shape == (10, 1)
        assert labels.shape == (10, 1)
    ########################################

    ########################################
    # lists with same lengths
    indices_list = [
        np.arange(0, adata.n_obs // 2),
        np.arange(adata.n_obs // 2, adata.n_obs),
    ]
    dataloader = ConcatDataLoader(
        manager, indices_list=indices_list, batch_size=10, drop_last=False
    )
    assert len(dataloader) == ceil((adata.n_obs / 2) / 10)
    dataloader = ConcatDataLoader(
        manager, indices_list=indices_list, batch_size=10, drop_last=True
    )
    assert len(dataloader) == floor((adata.n_obs / 2) / 10)

    for batch in dataloader:
        assert len(batch) == len(indices_list)
        X, batch, labels = (
            batch[0][REGISTRY_KEYS.X_KEY],
            batch[0][REGISTRY_KEYS.BATCH_KEY],
            batch[0][REGISTRY_KEYS.LABELS_KEY],
        )
        assert isinstance(X, torch.Tensor)
        assert isinstance(batch, torch.Tensor)
        assert isinstance(labels, torch.Tensor)
        assert X.shape == (10, n_genes)
        assert batch.shape == (10, 1)
        assert labels.shape == (10, 1)
    ########################################

    ########################################
    # lists with different lengths
    start1, end1 = 0, floor(adata.n_obs * (1 / 8))  # 1/8 of the data
    start2, end2 = end1, floor(adata.n_obs * (1 / 8 + 1 / 4))  # 1/4 of the data
    start3, end3 = end2, adata.n_obs  # 5/8 of the data
    batch_size = adata.n_obs // 8  # 1/8 of the data
    indices_list = [
        np.arange(start1, end1),  # dl1: [0]
        np.arange(start2, end2),  # dl2: [1][2]
        np.arange(start3, end3),  # dl3: [3][4][5][6][7]
    ]
    dataloader = ConcatDataLoader(
        manager, indices_list=indices_list, batch_size=batch_size, drop_last=False
    )
    assert len(dataloader) == ceil((adata.n_obs * (5 / 8)) / batch_size)
    dataloader = ConcatDataLoader(
        manager, indices_list=indices_list, batch_size=batch_size, drop_last=True
    )
    assert len(dataloader) == floor((adata.n_obs * (5 / 8)) / batch_size)

    # expected cycling
    # dl1: [0][0][0][0][0]
    # dl2: [1][2][1][2][1]
    # dl3: [3][4][5][6][7]
    prev_batch, prev_prev_batch = None, None
    for batch in dataloader:
        assert len(batch) == len(indices_list)
        for i, _batch in enumerate(batch):
            for key, value in _batch.items():
                assert isinstance(value, torch.Tensor)
                if i == 0 and prev_batch is not None:
                    assert torch.equal(value, prev_batch[i][key])
                elif i == 1 and prev_prev_batch is not None:
                    assert torch.equal(value, prev_prev_batch[i][key])

        prev_prev_batch = prev_batch
        prev_batch = batch
    ########################################
