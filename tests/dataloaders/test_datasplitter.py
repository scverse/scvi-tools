from __future__ import annotations

from math import ceil, floor

import numpy as np
import pytest
from sparse_utils import TestSparseModel

import scvi
from tests.data.utils import generic_setup_adata_manager


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


@pytest.mark.parametrize("sparse_format", ["csr_matrix", "csc_matrix"])
def test_datasplitter_load_sparse_tensor(
    sparse_format: str,
    accelerator: str,
    devices: list | str | int,
):
    adata = scvi.data.synthetic_iid(sparse_format=sparse_format)
    TestSparseModel.setup_anndata(adata)
    model = TestSparseModel(adata)
    model.train(
        accelerator=accelerator,
        devices=devices,
        expected_sparse_layout=sparse_format.split("_")[0],
    )
