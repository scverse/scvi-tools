from math import ceil, floor
from typing import Literal

import numpy as np
import pytest
import torch

import scvi
from tests.dataset.utils import generic_setup_adata_manager


class TestSparseDataSplitter(scvi.dataloaders.DataSplitter):
    def __init__(self, *args, expected_sparse_layout: Literal["csr", "csc"], **kwargs):
        if expected_sparse_layout == "csr":
            self.expected_sparse_layout = torch.sparse_csr
        elif expected_sparse_layout == "csc":
            self.expected_sparse_layout = torch.sparse_csc

        super().__init__(*args, **kwargs)

    def on_after_batch_transfer(self, batch, dataloader_idx):
        X = batch.get(scvi.REGISTRY_KEYS.X_KEY)
        assert isinstance(X, torch.Tensor)
        assert X.layout is self.expected_sparse_layout


def test_datasplitter_shuffle():
    adata = scvi.data.synthetic_iid()
    manager = generic_setup_adata_manager(adata)

    datasplitter = scvi.dataloaders.DataSplitter(
        manager, train_size=0.5, validation_size=0.3, shuffle_set_split=False
    )
    datasplitter.setup()
    assert isinstance(datasplitter.train_idx, np.ndarray)
    assert isinstance(datasplitter.val_idx, np.ndarray)
    assert isinstance(datasplitter.test_idx, np.ndarray)

    n_train = ceil(adata.n_obs * 0.5)
    n_val = floor(adata.n_obs * 0.3)
    n_test = adata.n_obs - n_train - n_val

    np.testing.assert_array_equal(
        datasplitter.train_idx,
        np.arange(n_val, n_val + n_train),
    )
    np.testing.assert_array_equal(
        datasplitter.val_idx,
        np.arange(n_val),
    )
    np.testing.assert_array_equal(
        datasplitter.test_idx,
        np.arange(n_val + n_train, n_val + n_train + n_test),
    )


@pytest.mark.cuda
@pytest.mark.parametrize(
    "sparse_format", ["csr_matrix", "csr_array", "csc_matrix", "csc_array"]
)
def test_datasplitter_load_sparse_tensor(sparse_format: str):
    adata = scvi.data.synthetic_iid(sparse_format=sparse_format)
    manager = generic_setup_adata_manager(adata)

    datasplitter = TestSparseDataSplitter(
        manager, expected_sparse_layout=sparse_format.split("_")[0]
    )
    datasplitter.setup()
    dataloader = datasplitter.train_dataloader()
    _ = next(iter(dataloader))
