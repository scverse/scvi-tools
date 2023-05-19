from math import ceil, floor

import numpy as np

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
