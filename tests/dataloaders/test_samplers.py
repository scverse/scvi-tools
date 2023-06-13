from math import ceil, floor

import numpy as np
import pytest

import scvi
from scvi.dataloaders import BatchDistributedSampler
from tests.dataset.utils import generic_setup_adata_manager


def test_batchdistributedsampler_init(
    batch_size: int = 128,
    n_batches: int = 2,
):
    adata = scvi.data.synthetic_iid(batch_size=batch_size, n_batches=n_batches)
    manager = generic_setup_adata_manager(adata)
    dataset = manager.create_torch_dataset()

    sampler = BatchDistributedSampler(
        dataset,
        num_replicas=1,
        rank=0,
        batch_size=batch_size,
        shuffle=True,
        drop_last=True,
        drop_dataset_tail=True,
    )
    assert sampler.batch_size == batch_size
    assert sampler.rank == 0
    assert sampler.num_replicas == 1
    assert sampler.shuffle
    assert sampler.drop_last_batch  # drop_last
    assert sampler.drop_last  # drop_dataset_tail


@pytest.mark.parametrize("drop_last", [True, False])
def test_batchdistributedsampler_drop_last(
    drop_last: bool,
    batch_size: int = 128,
    n_batches: int = 3,
    n_samplers: int = 2,
):
    adata = scvi.data.synthetic_iid(batch_size=batch_size, n_batches=n_batches)
    manager = generic_setup_adata_manager(adata)
    dataset = manager.create_torch_dataset()

    def check_samplers(samplers: list, sampler_batch_size: int):
        sampler_indices = [list(sampler) for sampler in samplers]
        n_obs_per_sampler = len(dataset) // n_samplers
        if drop_last:
            n_batches_per_sampler = floor(n_obs_per_sampler / sampler_batch_size)
        else:
            n_batches_per_sampler = ceil(n_obs_per_sampler / sampler_batch_size)

        assert len(sampler_indices) == n_samplers
        assert len(sampler_indices[0]) == len(sampler_indices[1])
        for indices in sampler_indices:
            assert len(indices) == n_batches_per_sampler
            if drop_last:
                assert all(
                    len(indices[i]) == sampler_batch_size
                    for i in range(n_batches_per_sampler)
                )
            else:
                assert all(
                    len(indices[i]) == sampler_batch_size
                    for i in range(n_batches_per_sampler - 1)
                )
                assert len(indices[-1]) == n_obs_per_sampler % sampler_batch_size

    samplers = [
        BatchDistributedSampler(
            dataset,
            num_replicas=n_samplers,
            rank=i,
            drop_last=drop_last,
            batch_size=batch_size,
        )
        for i in range(n_samplers)
    ]
    check_samplers(samplers, batch_size)

    samplers = [
        BatchDistributedSampler(
            dataset,
            num_replicas=n_samplers,
            rank=i,
            drop_last=drop_last,
            batch_size=batch_size - 1,
        )
        for i in range(n_samplers)
    ]
    check_samplers(samplers, batch_size - 1)

    samplers = [
        BatchDistributedSampler(
            dataset,
            num_replicas=n_samplers,
            rank=i,
            drop_last=drop_last,
            batch_size=batch_size + 1,
        )
        for i in range(n_samplers)
    ]
    check_samplers(samplers, batch_size + 1)


def test_batchdistributedsampler_non_overlapping(
    batch_size: int = 128,
    n_batches: int = 3,
    n_samplers: int = 2,
):
    adata = scvi.data.synthetic_iid(batch_size=batch_size, n_batches=n_batches)
    manager = generic_setup_adata_manager(adata)
    dataset = manager.create_torch_dataset()

    samplers = [
        BatchDistributedSampler(
            dataset,
            num_replicas=n_samplers,
            rank=i,
            batch_size=batch_size,
        )
        for i in range(n_samplers)
    ]
    sampler_indices = [list(sampler) for sampler in samplers]
    sampler_indices = [set(np.concatenate(indices)) for indices in sampler_indices]

    # check that indices are non-overlapping
    for i in range(n_samplers):
        for j in range(i + 1, n_samplers):
            assert len(sampler_indices[i].intersection(sampler_indices[j])) == 0
