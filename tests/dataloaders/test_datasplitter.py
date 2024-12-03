from __future__ import annotations

from math import ceil, floor

import numpy as np
import pytest
from sparse_utils import TestSparseModel
from tests.data.utils import generic_setup_adata_manager

import scvi


class TestDataSplitters:
    def test_datasplitter_shuffle(self):
        adata = scvi.data.synthetic_iid()
        manager = generic_setup_adata_manager(adata)

        with pytest.raises(ValueError) as excinfo:
            scvi.dataloaders.DataSplitter(
                manager,
                train_size=1.5,
                validation_size=0.3,
                shuffle_set_split=False,
            )
        assert str(excinfo.value) == "Invalid train_size. Must be: 0 < train_size <= 1"

        with pytest.raises(ValueError) as excinfo:
            scvi.dataloaders.DataSplitter(
                manager, train_size=0.5, validation_size=1.3, shuffle_set_split=False
            )
        assert str(excinfo.value) == "Invalid validation_size. Must be 0 <= validation_size < 1"

        with pytest.raises(ValueError) as excinfo:
            scvi.dataloaders.DataSplitter(
                manager, train_size=0.5, validation_size=0.8, shuffle_set_split=False
            )
        assert str(excinfo.value) == "train_size + validation_size must be between 0 and 1"

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

    def test_datasplitter_external(self):
        adata = scvi.data.synthetic_iid()
        manager = generic_setup_adata_manager(adata)

        # those can be inputs to the run
        valid_sec = 0.6
        test_sec = 0.8

        # random split of data
        train_ind, valid_ind, test_ind = np.split(
            [int(x) for x in adata.obs.sample(frac=1).index],
            [int(valid_sec * len(adata.obs)), int(test_sec * len(adata.obs))],
        )
        datasplitter = scvi.dataloaders.DataSplitter(
            manager, external_indexing=[train_ind, valid_ind, test_ind]
        )
        datasplitter.setup()
        assert isinstance(datasplitter.train_idx, np.ndarray)
        assert isinstance(datasplitter.val_idx, np.ndarray)
        assert isinstance(datasplitter.test_idx, np.ndarray)

        n_train = len(datasplitter.train_idx)
        n_val = len(datasplitter.val_idx)
        n_test = adata.n_obs - n_train - n_val

        # check for no ovarlapping and range
        assert n_test == len(datasplitter.test_idx)
        assert len(np.intersect1d(datasplitter.train_idx, datasplitter.val_idx)) == 0
        assert len(np.intersect1d(datasplitter.train_idx, datasplitter.test_idx)) == 0
        assert len(np.intersect1d(datasplitter.test_idx, datasplitter.val_idx)) == 0

    def test_datasplitter_external_with_overlap(self):
        adata = scvi.data.synthetic_iid()
        manager = generic_setup_adata_manager(adata)

        # those can be inputs to the run (this settings will also check for missing element)
        valid_sec = 0.6
        test_sec = 0.8

        # random split of data
        train_ind, valid_ind, test_ind = np.split(
            [int(x) for x in adata.obs.sample(frac=1).index],
            [int(valid_sec * len(adata.obs)), int(test_sec * len(adata.obs))],
        )
        valid_ind = np.append(valid_ind, test_ind[0])
        with pytest.raises(ValueError) as excinfo:
            scvi.dataloaders.DataSplitter(
                manager, external_indexing=[train_ind, valid_ind, test_ind]
            )
        assert str(excinfo.value) == "There are overlapping indexing between test and valid sets"

        train_ind = np.append(train_ind, test_ind[1])
        with pytest.raises(ValueError) as excinfo:
            scvi.dataloaders.DataSplitter(
                manager, external_indexing=[train_ind, valid_ind, test_ind]
            )
        assert str(excinfo.value) == "There are overlapping indexing between train and test sets"

        train_ind = np.append(train_ind, valid_ind[0])
        with pytest.raises(ValueError) as excinfo:
            scvi.dataloaders.DataSplitter(
                manager, external_indexing=[train_ind, valid_ind, test_ind]
            )
        assert str(excinfo.value) == "There are overlapping indexing between train and valid sets"

    def test_datasplitter_external_with_missing_indices(self):
        adata = scvi.data.synthetic_iid()
        manager = generic_setup_adata_manager(adata)

        # those can be inputs to the run (this settings will also check for missing element)
        valid_sec = 0.6
        test_sec = 1

        # random split of data
        train_ind, valid_ind, test_ind = np.split(
            [int(x) for x in adata.obs.sample(frac=1).index],
            [int(valid_sec * len(adata.obs)), int(test_sec * len(adata.obs))],
        )

        scvi.dataloaders.DataSplitter(manager, external_indexing=[train_ind, valid_ind, None])

        scvi.dataloaders.DataSplitter(manager, external_indexing=[train_ind, valid_ind])

        with pytest.raises(Warning) as excinfo:
            scvi.dataloaders.DataSplitter(manager, external_indexing=[train_ind])
        assert str(excinfo.value) == "There are missing indices please fix or remove those lines"

        with pytest.raises(Warning) as excinfo:
            scvi.dataloaders.DataSplitter(manager, external_indexing=[])
        assert str(excinfo.value) == "There are missing indices please fix or remove those lines"

        with pytest.raises(ValueError) as excinfo:
            scvi.dataloaders.DataSplitter(manager, external_indexing=train_ind)
        assert str(excinfo.value) == "External indexing is not of list type"

        with pytest.raises(ValueError) as excinfo:
            scvi.dataloaders.DataSplitter(manager, external_indexing=[[train_ind]])
        assert str(excinfo.value) == "One of the given external indexing arrays is not a np.array"

    def test_datasplitter_external_with_duplicates(self):
        adata = scvi.data.synthetic_iid()
        manager = generic_setup_adata_manager(adata)

        # those can be inputs to the run (this settings will also check for missing elements)
        valid_sec = 0.6
        test_sec = 0.8

        # random split of data
        train_ind, valid_ind, test_ind = np.split(
            [int(x) for x in adata.obs.sample(frac=1).index],
            [int(valid_sec * len(adata.obs)), int(test_sec * len(adata.obs))],
        )
        # Add duplicates
        test_ind = np.append(test_ind, test_ind[0])
        with pytest.raises(Warning) as excinfo:
            scvi.dataloaders.DataSplitter(
                manager, external_indexing=[train_ind, valid_ind, test_ind]
            )
        assert str(excinfo.value) == "There are duplicate indexing in test set"

        valid_ind = np.append(valid_ind, valid_ind[0])
        with pytest.raises(Warning) as excinfo:
            scvi.dataloaders.DataSplitter(manager, external_indexing=[train_ind, valid_ind])
        assert str(excinfo.value) == "There are duplicate indexing in valid set"

        train_ind = np.append(train_ind, train_ind[0])
        with pytest.raises(Warning) as excinfo:
            scvi.dataloaders.DataSplitter(manager, external_indexing=[train_ind])
        assert str(excinfo.value) == "There are duplicate indexing in train set"


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
