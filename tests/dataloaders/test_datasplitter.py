from __future__ import annotations

from math import ceil, floor

import numpy as np
import pytest
from sparse_utils import TestSparseModel
from tests.data.utils import generic_setup_adata_manager

import scvi


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

def test_datasplitter_external():
    adata = scvi.data.synthetic_iid()
    manager = generic_setup_adata_manager(adata)

    #those can be inputs to the run
    valid_sec = 0.6
    test_sec = 0.8

    #random split of data
    train_ind, validate_ind, test_ind = np.split([int(x) for x in adata.obs.sample(frac=1).index],
                                                 [int(valid_sec * len(adata.obs)),
                                                  int(test_sec * len(adata.obs))])
    datasplitter = scvi.dataloaders.DataSplitter(
        manager, use_external_indexing=True, external_indexing=[train_ind, validate_ind, test_ind]
    )
    datasplitter.setup()
    assert isinstance(datasplitter.train_idx, np.ndarray)
    assert isinstance(datasplitter.val_idx, np.ndarray)
    assert isinstance(datasplitter.test_idx, np.ndarray)

    n_train = len(datasplitter.train_idx)
    n_val = len(datasplitter.val_idx)
    n_test = adata.n_obs - n_train - n_val

    #check for no ovarlapping and range
    assert n_test==len(datasplitter.test_idx)
    assert len(np.intersect1d(datasplitter.train_idx, datasplitter.val_idx))==0
    assert len(np.intersect1d(datasplitter.train_idx, datasplitter.test_idx))==0
    assert len(np.intersect1d(datasplitter.test_idx, datasplitter.val_idx))==0


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
