from __future__ import annotations

from math import ceil, floor

import numpy as np
import pytest
from sparse_utils import TestSparseModel

import scvi
from tests.dataset.utils import generic_setup_adata_manager


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


@pytest.mark.parametrize(
    "sparse_format", ["csr_matrix", "csr_array", "csc_matrix", "csc_array"]
)
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


@pytest.mark.parametrize("shuffle_set_split", [False, True])
def test_contrastive_datasplitter(shuffle_set_split: bool):
    adata = scvi.data.synthetic_iid(n_batches=2)
    adata = adata[:-3, :]  # Unequal technical batch sizes.
    adata.layers["raw_counts"] = adata.X.copy()
    background_indices = (
        adata.obs.index[adata.obs["batch"] == "batch_0"].astype(int).tolist()
    )
    target_indices = (
        adata.obs.index[adata.obs["batch"] == "batch_1"].astype(int).tolist()
    )
    manager = generic_setup_adata_manager(
        adata=adata, batch_key="batch", labels_key="labels", layer="raw_counts"
    )
    contrastive_datasplitter = scvi.dataloaders.ContrastiveDataSplitter(
        adata_manager=manager,
        background_indices=background_indices,
        target_indices=target_indices,
        train_size=0.5,
        validation_size=0.3,
        shuffle_set_split=shuffle_set_split,
    )
    contrastive_datasplitter.setup()

    assert isinstance(contrastive_datasplitter.background_train_idx, list)
    assert isinstance(contrastive_datasplitter.background_val_idx, list)
    assert isinstance(contrastive_datasplitter.background_test_idx, list)
    assert isinstance(contrastive_datasplitter.target_train_idx, list)
    assert isinstance(contrastive_datasplitter.target_val_idx, list)
    assert isinstance(contrastive_datasplitter.target_test_idx, list)

    n_background = len(background_indices)
    n_background_train = ceil(n_background * 0.5)
    n_background_val = floor(n_background * 0.3)
    n_background_test = n_background - n_background_train - n_background_val

    assert len(contrastive_datasplitter.background_train_idx) == n_background_train
    assert len(contrastive_datasplitter.background_val_idx) == n_background_val
    assert len(contrastive_datasplitter.background_test_idx) == n_background_test

    n_target = len(target_indices)
    n_target_train = ceil(n_target * 0.5)
    n_target_val = floor(n_target * 0.3)
    n_target_test = n_target - n_target_train - n_target_val

    assert len(contrastive_datasplitter.target_train_idx) == n_target_train
    assert len(contrastive_datasplitter.target_val_idx) == n_target_val
    assert len(contrastive_datasplitter.target_test_idx) == n_target_test

    if not shuffle_set_split:
        np.testing.assert_array_equal(
            np.array(contrastive_datasplitter.background_val_idx),
            np.array(background_indices)[:n_background_val],
        )
        np.testing.assert_array_equal(
            np.array(contrastive_datasplitter.background_train_idx),
            np.array(background_indices)[
                n_background_val : (n_background_val + n_background_train)
            ],
        )
        np.testing.assert_array_equal(
            np.array(contrastive_datasplitter.background_test_idx),
            np.array(background_indices)[(n_background_val + n_background_train) :],
        )

        np.testing.assert_array_equal(
            np.array(contrastive_datasplitter.target_val_idx),
            np.array(target_indices)[:n_target_val],
        )
        np.testing.assert_array_equal(
            np.array(contrastive_datasplitter.target_train_idx),
            np.array(target_indices)[n_target_val : (n_target_val + n_target_train)],
        )
        np.testing.assert_array_equal(
            np.array(contrastive_datasplitter.target_test_idx),
            np.array(target_indices)[(n_target_val + n_target_train) :],
        )
