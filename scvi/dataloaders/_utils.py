import logging
from math import ceil, floor
from typing import List, Optional, Tuple

import numpy as np
from _docstrings import data_splitting_dsp

from scvi.utils._exceptions import InvalidParameterError

logger = logging.getLogger(__name__)


def _validate_data_split_sizes(
    n_obs: int,
    train_size: float,
    validation_size: Optional[float],
) -> Tuple[int, int, int]:
    if train_size > 1.0 or train_size <= 0.0:
        raise InvalidParameterError(
            "train_size",
            train_size,
            additional_message="`train_size` must be between 0 and 1.",
        )

    n_train = ceil(train_size * n_obs)

    if validation_size is None:
        n_validation = n_obs - n_train
    elif validation_size >= 1.0 or validation_size < 0.0:
        raise InvalidParameterError(
            "validation_size",
            validation_size,
            additional_message="`validation_size` must be between 0 and 1.",
        )
    elif (train_size + validation_size) > 1:
        raise InvalidParameterError(
            "train_size + validation_size",
            train_size + validation_size,
            additional_message="`train_size + validation_size` must be between 0 and 1.",
        )
    else:
        n_validation = floor(n_obs * validation_size)

    n_test = n_obs - n_train - n_validation

    logging.info(
        f"Using {n_train} observations for training, {n_validation} for validation "
        f"and {n_test} for testing."
    )
    return n_train, n_validation, n_test


def _validate_data_split_indices(
    n_obs: int,
    train_indices: List[int],
    validation_indices: Optional[List[int]],
) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    train_indices = np.array(train_indices)
    if validation_indices is not None:
        validation_indices = np.array(validation_indices)

    if np.amax(train_indices) >= n_obs or np.amin(train_indices) < 0:
        raise InvalidParameterError(
            "train_indices",
            train_indices,
            additional_message="`train_indices` contains invalid indices.",
        )

    if validation_indices is None:
        validation_indices = np.setdiff1d(np.arange(n_obs), train_indices)
    elif np.amax(validation_indices) >= n_obs or np.amin(validation_indices) < 0:
        raise InvalidParameterError(
            "validation_indices",
            validation_indices,
            additional_message="`validation_indices` contains invalid indices.",
        )

    union_indices = np.union1d(train_indices, validation_indices)
    test_indices = np.setdiff1d(np.arange(n_obs), union_indices)

    logging.info(
        f"Using {len(train_indices)} observations for training, "
        f"{len(validation_indices)} for validation and {len(test_indices)} for "
        "testing."
    )

    return train_indices, validation_indices, test_indices


@data_splitting_dsp.dedent
def validate_data_split(
    n_obs: int,
    train_size: Optional[float],
    validation_size: Optional[float],
    train_indices: Optional[List[int]],
    validation_indices: Optional[List[int]],
):
    """Validate data splitting parameters.

    Parameters
    ----------
    %(n_obs)s
    %(train_size)s
    %(validation_size)s
    %(train_indices)s
    %(validation_indices)s

    Returns
    -------
    * The number of training observations.
    * The number of validation observations.
    * The number of test observations.
    * The training observation indices if `train_indices` is not `None`.
    * The validation observation indices if `train_indices` is not `None`.te
    """
    if train_size is None and train_indices is None:
        raise ValueError("Either `train_size` or `train_indices` must be specified.")
    if train_size is not None and train_indices is not None:
        raise ValueError("`train_size` and `train_indices` cannot both be specified.")

    if train_size is None and validation_size is not None:
        raise ValueError(
            "`train_size` must be specified if `validation_size` is specified."
        )
    if train_indices is None and validation_indices is not None:
        raise ValueError(
            "`train_indices` must be specified if `validation_indices` is specified."
        )

    if train_size is not None:
        n_train, n_validation, n_test = _validate_data_split_sizes(
            n_obs, train_size, validation_size
        )
        return n_train, n_validation, n_test
    else:
        train_indices, validation_indices, test_indices = _validate_data_split_indices(
            n_obs, train_indices, validation_indices
        )
        return train_indices, validation_indices, test_indices
