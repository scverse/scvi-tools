import logging
from typing import TYPE_CHECKING, Dict, List, Union

import h5py
import numpy as np
import pandas as pd
from anndata._core.sparse_dataset import SparseDataset
from scipy.sparse import issparse
from torch.utils.data import Dataset

from scvi._constants import REGISTRY_KEYS

if TYPE_CHECKING:
    from ._manager import AnnDataManager

logger = logging.getLogger(__name__)


class AnnTorchDataset(Dataset):
    """Extension of torch dataset to get tensors from AnnData.

    Parameters
    ----------
    adata_manager
        :class:`~scvi.data.AnnDataManager` object with a registered AnnData object.
    getitem_tensors
        Dictionary with keys representing keys in data registry (``adata_manager.data_registry``)
        and value equal to desired numpy loading type (later made into torch tensor) or list of
        such keys. A list can be used to subset to certain keys in the event that more tensors than
        needed have been registered. If ``None``, defaults to all registered data.

    Examples
    --------
    >>> sd = AnnTorchDataset(adata_manager, ['X', 'batch'])
    >>> sd[[0, 1, 2]] # Return the X and batch data both as np.float32
    >>> sd = AnnTorchDataset(adata_manager, {'X':np.int64, 'batch':np.float32]})
    >>> sd[[0, 1, 2]] # Return the X as np.int64 and batch as np.float32
    """

    def __init__(
        self,
        adata_manager: "AnnDataManager",
        getitem_tensors: Union[List[str], Dict[str, type]] = None,
    ):
        if adata_manager.adata is None:
            raise ValueError(
                "Please run register_fields() on your AnnDataManager object first."
            )
        self.adata_manager = adata_manager
        self.is_backed = adata_manager.adata.isbacked
        self.attributes_and_types = None
        if getitem_tensors is not None:
            data_registry = adata_manager.data_registry
            for key in (
                getitem_tensors.keys()
                if isinstance(getitem_tensors, dict)
                else getitem_tensors
            ):
                if key not in data_registry.keys():
                    raise ValueError(
                        f"{key} required for model but not registered with AnnDataManager."
                    )
        self.getitem_tensors = getitem_tensors
        self._setup_getitem()
        self._set_data_attr()

    @property
    def registered_keys(self):
        """Returns the keys of the mappings in scvi data registry."""
        return self.adata_manager.data_registry.keys()

    def _set_data_attr(self):
        """Sets data attribute.

        Reduces number of times anndata needs to be accessed
        """
        self.data = {
            key: self.adata_manager.get_from_registry(key)
            for key, _ in self.attributes_and_types.items()
        }

    def _setup_getitem(self):
        """Sets up the __getitem__ function used by PyTorch.

        By default, getitem will return every single item registered in the scvi data registry
        and will attempt to infer the correct type. np.float32 for continuous values, otherwise np.int64.

        If you want to specify which specific tensors to return you can pass in a List of keys from
        the scvi data registry. If you want to speficy specific tensors to return as well as their
        associated types, then you can pass in a dictionary with their type.
        """
        registered_keys = self.registered_keys
        getitem_tensors = self.getitem_tensors
        if isinstance(getitem_tensors, List):
            keys = getitem_tensors
            keys_to_type = {key: np.float32 for key in keys}
        elif isinstance(getitem_tensors, Dict):
            keys = getitem_tensors.keys()
            keys_to_type = getitem_tensors
        elif getitem_tensors is None:
            keys = registered_keys
            keys_to_type = {key: np.float32 for key in keys}
        else:
            raise ValueError(
                "getitem_tensors invalid type. Expected: List[str] or Dict[str, type] or None"
            )
        for key in keys:
            if key not in registered_keys:
                raise KeyError(f"{key} not in data_registry")

        self.attributes_and_types = keys_to_type

    def __getitem__(self, idx: List[int]) -> Dict[str, np.ndarray]:
        """Get tensors in dictionary from anndata at idx."""
        data_numpy = {}

        if self.is_backed and not isinstance(idx, int):
            # need to sort idxs for h5py datasets
            idx = np.sort(idx)
        for key, dtype in self.attributes_and_types.items():
            cur_data = self.data[key]
            # for backed anndata
            if isinstance(cur_data, h5py.Dataset) or isinstance(
                cur_data, SparseDataset
            ):
                sliced_data = cur_data[idx]
                if issparse(sliced_data):
                    sliced_data = sliced_data.toarray()
                sliced_data = sliced_data.astype(dtype)
            elif isinstance(cur_data, np.ndarray):
                sliced_data = cur_data[idx].astype(dtype)
            elif isinstance(cur_data, pd.DataFrame):
                sliced_data = cur_data.iloc[idx, :].to_numpy().astype(dtype)
            elif issparse(cur_data):
                sliced_data = cur_data[idx].toarray().astype(dtype)
            # for minified  anndata, we need this because we can have a string
            # cur_data, which is the value of the MINIFY_TYPE_KEY in adata.uns,
            # used to record the type data minification
            # TODO: Adata manager should have a list of which fields it will load
            elif isinstance(cur_data, str) and key == REGISTRY_KEYS.MINIFY_TYPE_KEY:
                continue
            else:
                raise TypeError(f"{key} is not a supported type")
            # Make a row vector if only one element is selected
            # this is because our dataloader disables automatic batching
            # Normally, this would be handled by the default collate fn
            if isinstance(idx, int) or len(idx) == 1:
                sliced_data = sliced_data.reshape(1, -1)
            data_numpy[key] = sliced_data

        return data_numpy

    def get_data(self, scvi_data_key: str):
        """Get the tensor associated with a key."""
        tensors = self.__getitem__(idx=list(range(self.__len__())))
        return tensors[scvi_data_key]

    def __len__(self):
        return self.adata_manager.adata.shape[0]
