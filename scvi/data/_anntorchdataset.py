from __future__ import annotations

import logging
from typing import TYPE_CHECKING

import h5py
import numpy as np
import pandas as pd
import torch

try:
    from anndata._core.sparse_dataset import SparseDataset
except ImportError:
    # anndata >= 0.10.0
    from anndata._core.sparse_dataset import (
        BaseCompressedSparseDataset as SparseDataset,
    )

from scipy.sparse import issparse
from torch.utils.data import Dataset

from scvi._constants import REGISTRY_KEYS
from scvi.utils._exceptions import InvalidParameterError

if TYPE_CHECKING:
    from ._manager import AnnDataManager
from ._utils import registry_key_to_default_dtype, scipy_to_torch_sparse

logger = logging.getLogger(__name__)


class AnnTorchDataset(Dataset):
    """Extension of :class:`~torch.utils.data.Dataset` for :class:`~anndata.AnnData` objects.

    Parameters
    ----------
    adata_manager
        :class:`~scvi.data.AnnDataManager` object with a registered AnnData object.
    getitem_tensors
        Specifies the keys in the data registry (``adata_manager.data_registry``) to return in
        ``__getitem__``. One of the following:

        * ``dict``: Keys correspond to keys in the data registry and values correspond to the
        desired :class:`~np.dtype` of the returned data.
        * ``list``: Elements correspond to keys in the data registry. Continuous data will be
        returned as :class:`~np.float32` and discrete data will be returned as :class:`~np.int64`.
        * ``None``: All registered data will be returned. Continuous data will be returned as
        :class:`~np.float32` and discrete data will be returned as :class:`~np.int64`.
    load_sparse_tensor
        `EXPERIMENTAL` If ``True``, loads data with sparse CSR or CSC layout as a
        :class:`~torch.Tensor` with the same layout. Can lead to speedups in data transfers to GPUs,
        depending on the sparsity of the data.
    """

    def __init__(
        self,
        adata_manager: AnnDataManager,
        getitem_tensors: list | dict[str, type] | None = None,
        load_sparse_tensor: bool = False,
    ):
        super().__init__()

        if adata_manager.adata is None:
            raise ValueError(
                "Please run ``register_fields`` on ``adata_manager`` first."
            )
        self.adata_manager = adata_manager
        self.keys_and_dtypes = getitem_tensors
        self.load_sparse_tensor = load_sparse_tensor

    @property
    def registered_keys(self):
        """Keys in the data registry."""
        return self.adata_manager.data_registry.keys()

    @property
    def keys_and_dtypes(self):
        """Keys and corresponding :class:`~np.dtype` of data to fetch in ``__getitem__``."""
        return self._keys_and_dtypes

    @keys_and_dtypes.setter
    def keys_and_dtypes(self, getitem_tensors: list | dict[str, type] | None):
        """Set keys and corresponding :class:`~np.dtype` of data to fetch in ``__getitem__``.

        Raises an error if any of the keys are not in the data registry.
        """
        if isinstance(getitem_tensors, list):
            keys_to_dtypes = {
                key: registry_key_to_default_dtype(key) for key in getitem_tensors
            }
        elif isinstance(getitem_tensors, dict):
            keys_to_dtypes = getitem_tensors
        elif getitem_tensors is None:
            keys_to_dtypes = {
                key: registry_key_to_default_dtype(key) for key in self.registered_keys
            }
        else:
            raise InvalidParameterError(
                param="getitem_tensors",
                value=getitem_tensors.__class__,
                valid=[list, dict, None],
            )

        for key in keys_to_dtypes:
            if key not in self.registered_keys:
                raise KeyError(f"{key} not found in the data registry.")

        self._keys_and_dtypes = keys_to_dtypes

    @property
    def data(self):
        """Dictionary of data tensors.

        First time this is accessed, data is fetched from the underlying
        :class:`~anndata.AnnData` object. Subsequent accesses will return the
        cached dictionary.
        """
        if not hasattr(self, "_data"):
            self._data = {
                key: self.adata_manager.get_from_registry(key)
                for key in self.keys_and_dtypes
            }
        return self._data

    def __len__(self):
        return self.adata_manager.adata.shape[0]

    def __getitem__(
        self, indexes: int | list[int] | slice
    ) -> dict[str, np.ndarray | torch.Tensor]:
        """Fetch data from the :class:`~anndata.AnnData` object.

        Parameters
        ----------
        indexes
            Indexes of the observations to fetch. Can be a single index, a list of indexes, or a
            slice.

        Returns
        -------
        Mapping of data registry keys to arrays of shape ``(n_obs, ...)``.
        """
        if isinstance(indexes, int):
            indexes = [indexes]  # force batched single observations

        if self.adata_manager.adata.isbacked and isinstance(indexes, (list, np.ndarray)):
            # need to sort indexes for h5py datasets
            indexes = np.sort(indexes)

        data_map = {}

        for key, dtype in self.keys_and_dtypes.items():
            data = self.data[key]

            if isinstance(data, (np.ndarray, h5py.Dataset)):
                sliced_data = data[indexes].astype(dtype, copy=False)
            elif isinstance(data, pd.DataFrame):
                sliced_data = data.iloc[indexes, :].to_numpy().astype(dtype, copy=False)
            elif issparse(data) or isinstance(data, SparseDataset):
                sliced_data = data[indexes].astype(dtype, copy=False)
                if self.load_sparse_tensor:
                    sliced_data = scipy_to_torch_sparse(sliced_data)
                else:
                    sliced_data = sliced_data.toarray()
            elif isinstance(data, str) and key == REGISTRY_KEYS.MINIFY_TYPE_KEY:
                # for minified  anndata, we need this because we can have a string
                # for `data``, which is the value of the MINIFY_TYPE_KEY in adata.uns,
                # used to record the type data minification
                # TODO: Adata manager should have a list of which fields it will load
                continue
            else:
                raise TypeError(f"{key} is not a supported type")

            data_map[key] = sliced_data

        return data_map
