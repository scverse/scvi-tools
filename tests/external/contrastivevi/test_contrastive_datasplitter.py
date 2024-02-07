from __future__ import annotations

from math import ceil, floor

import numpy as np
import pytest

import scvi


@pytest.mark.parametrize("shuffle_set_split", [False, True])
def test_contrastive_datasplitter(
    shuffle_set_split: bool,
    mock_contrastive_adata_manager,
    mock_background_indices,
    mock_target_indices,
):
    background_indices = mock_background_indices
    target_indices = mock_target_indices

    contrastive_datasplitter = scvi.external.contrastivevi.ContrastiveDataSplitter(
        adata_manager=mock_contrastive_adata_manager,
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
