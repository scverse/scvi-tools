import logging
from dataclasses import dataclass
from math import ceil, floor
from typing import List, Optional, Tuple, Union

import numpy as np
from anndata._core.sparse_dataset import SparseDataset
from h5py import Dataset
from pandas import DataFrame
from scipy.sparse import issparse
from torch import Tensor

from scvi import settings
from scvi.utils._exceptions import InvalidParameterError

from ._docstrings import datasplitter_dsp

logger = logging.getLogger(__name__)
ArrayLike = Union[np.ndarray, DataFrame, Dataset, SparseDataset, Tensor]


@dataclass
class SplitIndices:
    """Dataclass for storing the indices of the data split.

    Enforces that the indices are unique and supports addition.

    Parameters
    ----------
    train
        Indices of the training set.
    validation
        Indices of the validation set.
    test
        Indices of the test set.
    """

    train: List[int]
    validation: List[int]
    test: List[int]

    def __post_init__(self):
        self.train = np.array(self.train, dtype=int)
        self.validation = np.array(self.validation, dtype=int)
        self.test = np.array(self.test, dtype=int)

    def __add__(self, other):
        return SplitIndices(
            train=np.union1d(self.train, other.train),
            validation=np.union1d(self.validation, other.validation),
            test=np.union1d(self.test, other.test),
        )


def _validate_data_split_sizes(
    *,
    all_indices: np.ndarray,
    train_size: float,
    validation_size: Optional[float],
) -> Tuple[int, int, int]:
    if train_size > 1.0 or train_size <= 0.0:
        raise InvalidParameterError(
            "train_size",
            train_size,
            additional_message="`train_size` must be between 0 and 1.",
        )

    n_obs = len(all_indices)
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

    return n_train, n_validation, n_obs - (n_train + n_validation)


def _validate_data_split_indices(
    *,
    all_indices: np.ndarray,
    train_indices: List[int],
    validation_indices: Optional[List[int]],
) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    train_indices = np.unique(train_indices)
    if np.intersect1d(all_indices, train_indices).size != len(train_indices):
        raise InvalidParameterError(
            "train_indices",
            train_indices,
            additional_message="`train_indices` contains invalid indices.",
        )
    if validation_indices is None:
        validation_indices = np.setdiff1d(all_indices, train_indices)
    else:
        validation_indices = np.unique(validation_indices)
        if np.intersect1d(all_indices, validation_indices).size != len(
            validation_indices
        ):
            raise InvalidParameterError(
                "validation_indices",
                validation_indices,
                additional_message="`validation_indices` contains invalid indices.",
            )

    union_indices = np.union1d(train_indices, validation_indices)
    test_indices = np.setdiff1d(all_indices, union_indices)

    return train_indices, validation_indices, test_indices


def _make_data_split(
    *,
    all_indices: np.ndarray,
    n_train: int,
    n_validation: int,
    n_test: int,
    shuffle: bool,
) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    if shuffle:
        random_state = np.random.default_rng(seed=settings.seed)
        all_indices = random_state.permutation(all_indices)

    n_val_train = n_train + n_validation
    train_indices = all_indices[:n_train]
    validation_indices = all_indices[n_train:n_val_train]
    test_indices = all_indices[n_val_train:]

    return train_indices, validation_indices, test_indices


@datasplitter_dsp.dedent
def validate_data_split(
    *,
    n_obs: Optional[int] = None,
    all_indices: Optional[List[int]] = None,
    train_size: Optional[float] = None,
    validation_size: Optional[float] = None,
    train_indices: Optional[List[int]] = None,
    validation_indices: Optional[List[int]] = None,
    shuffle,
) -> SplitIndices:
    """Validate data splitting parameters.

    Parameters
    ----------
    %(param_n_obs)s
    %(param_all_indices)s
    %(param_train_size)s
    %(param_validation_size)s
    %(param_train_indices)s
    %(param_validation_indices)s

    Returns
    -------
    * The training observation indices.
    * The validation observation indices.
    * The test observation indices.
    """
    if n_obs is None and all_indices is None:
        raise ValueError("Either `n_obs` or `all_indices` must be specified.")
    if n_obs is not None and all_indices is not None:
        raise ValueError("`n_obs` and `all_indices` cannot both be specified.")

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

    if all_indices is None:
        all_indices = np.arange(n_obs)

    if train_size is not None:
        n_train, n_validation, n_test = _validate_data_split_sizes(
            all_indices=all_indices,
            train_size=train_size,
            validation_size=validation_size,
        )
        train_indices, validation_indices, test_indices = _make_data_split(
            all_indices=all_indices,
            n_train=n_train,
            n_validation=n_validation,
            n_test=n_test,
            shuffle=shuffle,
        )
    else:
        train_indices, validation_indices, test_indices = _validate_data_split_indices(
            all_indices=all_indices,
            train_indices=train_indices,
            validation_indices=validation_indices,
        )

    logging.info(
        f"Using {len(train_indices)} observations for training, "
        f"{len(validation_indices)} for validation and {len(test_indices)} for "
        "testing."
    )

    return SplitIndices(
        train=train_indices, validation=validation_indices, test=test_indices
    )


def slice_and_convert(
    data: ArrayLike,
    indices: Optional[List[int]] = None,
    dtype: Optional[str] = None,
) -> np.ndarray:
    """Slices and converts to the specified numpy dtype."""
    indices = indices or np.arange(len(data))

    if isinstance(data, Dataset) or isinstance(data, SparseDataset):
        _data = data[indices]
        if issparse(_data):
            _data = _data.toarray()
    elif isinstance(data, np.ndarray):
        _data = data[indices]
    elif isinstance(data, DataFrame):
        _data = data.iloc[indices, :].to_numpy()
    elif issparse(data):
        _data = data[indices].toarray()
    else:
        raise InvalidParameterError(
            param="data", value=data.__class__.__name__, valid=ArrayLike
        )

    if dtype is None:
        return _data
    return _data.astype(dtype)
