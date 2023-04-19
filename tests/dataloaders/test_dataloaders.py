from itertools import product
from math import ceil, floor

import numpy as np
import pytest
import torch
from torch.utils.data import SequentialSampler

from scvi import REGISTRY_KEYS
from scvi.data import synthetic_iid
from scvi.dataloaders import (
    AnnDataLoader,
    AnnTorchDataset,
    ConcatDataLoader,
    SemiSupervisedDataLoader,
)
from tests.dataset.utils import generic_setup_adata_manager, scanvi_setup_adata_manager


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


# TODO: Add disk backed, mudata tests


@pytest.mark.parametrize(
    "batch_size, n_batches, n_genes, sparse",
    product([128], [1, 2], [50], [True, False]),
)
def test_basic_concatdataloader(
    batch_size: int, n_batches: int, n_genes: int, sparse: bool
):
    adata = synthetic_iid(
        batch_size=batch_size,
        n_batches=n_batches,
        n_genes=n_genes,
        sparse=sparse,
    )
    manager = generic_setup_adata_manager(adata, batch_key="batch", labels_key="labels")
    n_obs = adata.n_obs

    # one list of indices, works the same as AnnDataLoader except needing to index batch
    indices_list = [np.arange(0, n_obs // 2)]
    dl = ConcatDataLoader(manager, indices_list, batch_size=10, drop_last=True)
    assert len(dl) == floor((n_obs / 2) / 10)  # floor since drop_last=True
    dl = ConcatDataLoader(manager, indices_list, batch_size=10, drop_last=False)
    assert len(dl) == ceil((n_obs / 2) / 10)  # ceil since drop_last=False

    batch = next(iter(dl))
    assert len(batch) == len(indices_list)
    for key, value in batch[0].items():
        assert key in (
            REGISTRY_KEYS.X_KEY,
            REGISTRY_KEYS.BATCH_KEY,
            REGISTRY_KEYS.LABELS_KEY,
        )
        assert isinstance(value, torch.Tensor)
        assert value.shape[0] == 10
        assert value.shape[1] == n_genes if key == REGISTRY_KEYS.X_KEY else 1

    # multiple lists of indices with same lengths
    indices_list = [np.arange(0, n_obs // 2), np.arange(n_obs // 2, n_obs)]
    dl = ConcatDataLoader(manager, indices_list, batch_size=10, drop_last=True)
    assert len(dl) == floor((n_obs / 2) / 10)  # floor since drop_last=True
    dl = ConcatDataLoader(manager, indices_list, batch_size=10, drop_last=False)
    assert len(dl) == ceil((n_obs / 2) / 10)  # ceil since drop_last=False

    batch = next(iter(dl))
    assert len(batch) == len(indices_list)
    for _batch in batch:
        for key, value in _batch.items():
            assert key in (
                REGISTRY_KEYS.X_KEY,
                REGISTRY_KEYS.BATCH_KEY,
                REGISTRY_KEYS.LABELS_KEY,
            )
            assert isinstance(value, torch.Tensor)
            assert value.shape[0] == 10
            assert value.shape[1] == n_genes if key == REGISTRY_KEYS.X_KEY else 1

    # multiple lists of indices with different lengths
    start1, end1 = 0, floor(n_obs * (1 / 8))  # 1/8 of the data
    start2, end2 = end1, floor(n_obs * (1 / 8 + 1 / 4))  # 1/4 of the data
    start3, end3 = end2, n_obs  # 5/8 of the data
    batch_size = n_obs // 8  # 1/8 of the data
    indices_list = [
        np.arange(start1, end1),  # dl1: [0]
        np.arange(start2, end2),  # dl2: [1][2]
        np.arange(start3, end3),  # dl3: [3][4][5][6][7]
    ]

    dl = ConcatDataLoader(manager, indices_list, batch_size=batch_size, drop_last=True)
    assert len(dl) == floor(
        (n_obs * (5 / 8)) / batch_size
    )  # floor since drop_last=True
    dl = ConcatDataLoader(manager, indices_list, batch_size=batch_size, drop_last=False)
    assert len(dl) == ceil((n_obs * (5 / 8)) / batch_size)  # ceil since drop_last=False

    # test that smaller dataloaders cycle correctly, expected cycling:
    # dl1: [0][0][0][0][0]
    # dl2: [1][2][1][2][1]
    # dl3: [3][4][5][6][7]
    prev_batch, prev_prev_batch = None, None
    for batch in dl:
        assert len(batch) == len(indices_list)
        for i, _batch in enumerate(batch):
            for key, value in _batch.items():
                assert key in (
                    REGISTRY_KEYS.X_KEY,
                    REGISTRY_KEYS.BATCH_KEY,
                    REGISTRY_KEYS.LABELS_KEY,
                )
                assert isinstance(value, torch.Tensor)
                assert value.shape[0] == batch_size
                assert value.shape[1] == n_genes if key == REGISTRY_KEYS.X_KEY else 1

                if i == 0 and prev_batch is not None:
                    assert torch.equal(value, prev_batch[i][key])
                elif i == 1 and prev_prev_batch is not None:
                    assert torch.equal(value, prev_prev_batch[i][key])

        prev_prev_batch = prev_batch
        prev_batch = batch


@pytest.mark.parametrize(
    "batch_size, n_batches, n_genes, n_labels, sparse",
    product([128], [1, 2], [50], [0, 1, 2, 3], [True, False]),
)
def test_basic_semisuperviseddataloader(
    batch_size: int, n_batches: int, n_genes: int, n_labels: int, sparse: bool
):
    adata = synthetic_iid(
        batch_size=batch_size,
        n_batches=n_batches,
        n_genes=n_genes,
        n_labels=n_labels,
        sparse=sparse,
    )
    labels_key = "labels" if n_labels > 0 else None
    unlabeled_category = "label_0" if n_labels > 2 else None
    label_counts = None
    if labels_key is not None:
        label_counts = adata.obs[labels_key].value_counts().to_dict()
    manager = scanvi_setup_adata_manager(
        adata, unlabeled_category=unlabeled_category, labels_key=labels_key
    )
    n_obs = adata.n_obs
    n_labeled = 0
    if labels_key is not None:
        n_labeled = sum([v for k, v in label_counts.items() if k != unlabeled_category])

    dl = SemiSupervisedDataLoader(manager, batch_size=128, drop_last=True)
    assert len(dl) == floor(n_obs / 128)  # floor since drop_last=True
    if labels_key is not None and n_labeled < 128:
        with pytest.raises(StopIteration):
            _ = next(iter(dl))
    else:
        batch = next(iter(dl))
        assert len(batch) == 2

    dl = SemiSupervisedDataLoader(manager, batch_size=128, drop_last=False)
    assert len(dl) == ceil(n_obs / 128)  # ceil since drop_last=False
    batch = next(iter(dl))
    _, labeled_obs = batch
    if labels_key is not None and n_labeled < 128:
        print(labeled_obs)
        assert labeled_obs[REGISTRY_KEYS.X_KEY].shape[0] == n_labeled


@pytest.mark.parametrize(
    "batch_size, n_batches, n_genes, n_labels, sparse",
    product([128], [1, 2], [50], [0, 1, 2, 3], [True, False]),
)
def test_n_samples_per_label_semisuperviseddataloader(
    batch_size: int, n_batches: int, n_genes: int, n_labels: int, sparse: bool
):
    # TODO
    pass
