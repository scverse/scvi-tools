from itertools import product

import pytest

from scvi.data import synthetic_iid
from scvi.dataloaders import AnnDataLoader
from tests.dataset.utils import generic_setup_adata_manager


@pytest.mark.parametrize(
    "batch_size, n_batches, sparse", product([128, 256], [1, 2], [True, False])
)
def test_basic_anndataloader(batch_size: int, n_batches: int, sparse: bool):
    adata = synthetic_iid(batch_size=batch_size, n_batches=n_batches, sparse=sparse)
    manager = generic_setup_adata_manager(adata, batch_key="batch", labels_key="labels")

    dl = AnnDataLoader(manager, batch_size=batch_size)
    assert len(dl) == n_batches
    batch = next(iter(dl))
    print(batch)


@pytest.mark.cuda
@pytest.mark.parametrize(
    "batch_size, n_batches, sparse", product([128, 256], [1, 2], [True, False])
)
def test_cuda_anndataloader(cuda: bool, batch_size: int, n_batches: int, sparse: bool):
    synthetic_iid(batch_size=batch_size, n_batches=n_batches, sparse=sparse)


def test_basic_concatdataloader():
    pass


def test_basic_semisuperviseddataloader():
    pass
