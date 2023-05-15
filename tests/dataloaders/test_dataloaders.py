import pytest
from torch.utils.data import BatchSampler, RandomSampler, SequentialSampler

from scvi.data import synthetic_iid
from scvi.dataloaders import (
    AnnDataLoader,
    AnnTorchDataset,
)
from tests.dataset.utils import generic_setup_adata_manager


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
