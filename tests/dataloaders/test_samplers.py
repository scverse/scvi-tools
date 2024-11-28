from math import ceil, floor

import numpy as np
import os
import pytest
from tests.data.utils import generic_setup_adata_manager

import scvi
from scvi.dataloaders import BatchDistributedSampler


@pytest.mark.parametrize("num_processes", [1, 2])
def test_batchdistributedsampler_init(
    num_processes: int,
    save_path: str,
    batch_size: int = 128,
    n_batches: int = 2,
):
    adata = scvi.data.synthetic_iid(batch_size=batch_size, n_batches=n_batches)
    manager = generic_setup_adata_manager(adata)
    dataset = manager.create_torch_dataset()

    file_path = save_path + "/dist_file"
    if os.path.exists(file_path):  # Check if the file exists
        os.remove(file_path)

    sampler = BatchDistributedSampler(
        dataset,
        num_replicas=1,
        rank=0,
        batch_size=batch_size,
        shuffle=True,
        drop_last=True,
        drop_dataset_tail=True,
        num_processes=num_processes,
        save_path=save_path
    )
    assert sampler.batch_size == batch_size
    assert sampler.rank == 0
    assert sampler.num_replicas == 1
    assert sampler.shuffle
    assert sampler.drop_last_batch  # drop_last
    assert sampler.drop_last  # drop_dataset_tail


@pytest.mark.parametrize("drop_last", [True, False])
@pytest.mark.parametrize("drop_dataset_tail", [True, False])
@pytest.mark.parametrize("num_processes", [1, 2])
def test_batchdistributedsampler_drop_last(
    num_processes: int,
    drop_last: bool,
    drop_dataset_tail: bool,
    save_path: str,
    batch_size: int = 128,
    n_batches: int = 3,
    num_replicas: int = 2,
):
    """Expected behavior:

    drop_last: False
    drop_dataset_tail: False
    - samplers see m = ceil(n_obs / num_replicas) observations
    - samplers have b = ceil(m / batch_size) batches

    drop_last: True
    drop_dataset_tail: False
    - samplers see m = ceil(n_obs / num_replicas) observations
    - samplers have b = floor(m / batch_size) batches
    - error if b == 0

    drop_last: False
    drop_dataset_tail: True
    - samplers see m = floor(n_obs / num_replicas) observations
    - samplers have b = ceil(m / batch_size) batches
    - error if m == 0

    drop_last: True
    drop_dataset_tail: True
    - samplers see m = floor(n_obs / num_replicas) observations
    - samplers have b = floor(m / batch_size) batches
    - error if m == 0 or b == 0
    """
    adata = scvi.data.synthetic_iid(batch_size=batch_size, n_batches=n_batches)
    manager = generic_setup_adata_manager(adata)
    dataset = manager.create_torch_dataset()

    def check_samplers(samplers: list, sampler_batch_size: int):
        sampler_indices = [list(sampler) for sampler in samplers]

        if drop_dataset_tail:
            n_obs_per_sampler = floor(len(dataset) / num_replicas)
        else:
            n_obs_per_sampler = ceil(len(dataset) / num_replicas)

        if drop_last:
            n_batches_per_sampler = floor(n_obs_per_sampler / sampler_batch_size)
            batch_sizes = [sampler_batch_size] * n_batches_per_sampler
        else:
            n_batches_per_sampler = ceil(n_obs_per_sampler / sampler_batch_size)
            if n_obs_per_sampler % sampler_batch_size == 0:
                batch_sizes = [sampler_batch_size] * n_batches_per_sampler
            else:
                batch_sizes = [sampler_batch_size] * (n_batches_per_sampler - 1) + [
                    n_obs_per_sampler % sampler_batch_size
                ]

        effective_n_obs_per_sampler = sum(batch_sizes)

        assert len(sampler_indices) == num_replicas
        for batch_indices in sampler_indices:
            all_indices = np.concatenate(batch_indices)

            assert len(batch_indices) == n_batches_per_sampler
            assert len(all_indices) == effective_n_obs_per_sampler
            assert [len(indices) for indices in batch_indices] == batch_sizes

    file_path = save_path + "/dist_file"
    if os.path.exists(file_path):  # Check if the file exists
        os.remove(file_path)

    for sampler_batch_size in [batch_size, batch_size - 1, batch_size + 1]:
        samplers = [
            BatchDistributedSampler(
                dataset,
                num_replicas=num_replicas,
                rank=i,
                batch_size=sampler_batch_size,
                drop_last=drop_last,
                drop_dataset_tail=drop_dataset_tail,
                num_processes=num_processes,
                save_path=save_path
            )
            for i in range(num_replicas)
        ]
        check_samplers(samplers, sampler_batch_size)


@pytest.mark.parametrize("num_processes", [1, 2])
def test_batchdistributedsampler_indices(
    num_processes: int,
    save_path: str,
    batch_size: int = 128,
    n_batches: int = 3,
    num_replicas: int = 2,
):
    adata = scvi.data.synthetic_iid(batch_size=batch_size, n_batches=n_batches)
    manager = generic_setup_adata_manager(adata)
    dataset = manager.create_torch_dataset()

    file_path = save_path + "/dist_file"
    if os.path.exists(file_path):  # Check if the file exists
        os.remove(file_path)

    samplers = [
        BatchDistributedSampler(
            dataset,
            num_replicas=num_replicas,
            rank=i,
            batch_size=batch_size,
            num_processes=num_processes,
            save_path=save_path
        )
        for i in range(num_replicas)
    ]
    sampler_indices = [list(sampler) for sampler in samplers]
    sampler_indices = [set(np.concatenate(indices)) for indices in sampler_indices]

    # check that indices are non-overlapping
    for i in range(num_replicas):
        for j in range(i + 1, num_replicas):
            assert len(sampler_indices[i].intersection(sampler_indices[j])) == 0

    # check that all indices are covered
    covered_indices = np.concatenate([np.array(list(indices)) for indices in sampler_indices])
    assert len(covered_indices) == len(dataset)
