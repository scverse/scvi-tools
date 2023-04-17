import logging
from typing import Dict, List, Optional, Union

import numpy as np
import torch
from torch.utils.data import Dataset

from scvi._constants import REGISTRY_KEYS
from scvi.data import AnnDataManager
from scvi.model._utils import parse_device_args

from ._docstrings import dataset_dsp
from ._utils import slice_and_convert

logger = logging.getLogger(__name__)


@dataset_dsp.dedent
class AnnTorchDataset(Dataset):
    """%(summary)s

    Parameters
    ----------
    %(param_adata_manager)s
    %(param_getitem_tensors)s
    %(param_accelerator)s
    %(param_device)s
    %(param_device_backed)s

    Examples
    --------
    >>> sd = AnnTorchDataset(adata_manager, ['X', 'batch'])
    >>> sd[[0, 1, 2]] # Return the X and batch data both as np.float32
    >>> sd = AnnTorchDataset(adata_manager, {'X':np.int64, 'batch':np.float32]})
    >>> sd[[0, 1, 2]] # Return the X as np.int64 and batch as np.float32
    """

    def __init__(
        self,
        adata_manager: AnnDataManager,
        getitem_tensors: Optional[Union[List[str], Dict[str, type]]] = None,
        accelerator: str = "auto",
        device: Union[int, str] = "auto",
        device_backed: bool = False,
    ):
        if adata_manager.adata is None:
            raise ValueError(
                "Please run `register_fields` on your `AnnDataManager` first."
            )
        self.adata_manager = adata_manager
        self.attributes_and_types = getitem_tensors

        _, _, self.device = parse_device_args(
            accelerator=accelerator,
            devices=device,
            return_device="torch",
            validate_single_device=True,
        )
        self.device_backed = device_backed
        self.backed_adata = adata_manager.adata.isbacked
        _ = self.data  # loads data from adata

    @property
    def registered_keys(self):
        """Returns the keys of the mappings in scvi data registry."""
        return self.adata_manager.data_registry.keys()

    @property
    def attributes_and_types(self) -> dict:
        """Returns the attributes and types to be loaded by `__getitem__`."""
        return self._attributes_and_types

    @attributes_and_types.setter
    def attributes_and_types(self, value: Optional[Union[list, dict]]):
        """Sets the attributes and types to be loaded by `__getitem__`."""
        if not isinstance(value, (list, dict, type(None))):
            raise ValueError(
                "`getitem_tensors` invalid type, expected a `list`, `dict`, or `None`."
            )
        if value is not None:
            for key in value:
                if key in self.registered_keys:
                    continue
                raise ValueError(
                    f"{key} required for the model but not registered with "
                    "`AnnDataManager`."
                )
        else:
            value = self.registered_keys

        if isinstance(value, dict):
            for key, dtype in value.items():
                if not isinstance(dtype, type):
                    raise ValueError(
                        f"`{key}` must have a valid `type` as its value, not {dtype}."
                    )

        if isinstance(value, list):
            value = {key: np.float32 for key in value}

        self._attributes_and_types = value

    @property
    def data(self) -> dict:
        """Sets the data attribute of the dataset."""
        if hasattr(self, "_data"):
            return self._data

        data = {}
        for key, dtype in self.attributes_and_types.items():
            _data = self.adata_manager.get_from_registry(key)

            if self.device_backed and key != REGISTRY_KEYS.MINIFY_TYPE_KEY:
                _data = slice_and_convert(_data, dtype=dtype)
                _data = torch.from_numpy(_data).to(self.device)

            data[key] = _data

        self._data = data
        return self._data

    def __getitem__(
        self, indices: Union[list, int]
    ) -> Dict[str, Union[np.ndarray, torch.Tensor]]:
        """Slice data attributes at the specified indices."""
        if isinstance(indices, slice):
            indices = np.arange(*indices.indices(len(self)))
        elif isinstance(indices, int) or isinstance(indices, np.integer):
            indices = np.array([indices])
        elif isinstance(indices, np.ndarray):
            indices = indices.astype(int)
        elif isinstance(indices, list):
            indices = np.array(indices)
        else:
            raise TypeError(
                "`indices` must be a `slice`, `int`, `np.ndarray`, or `list`."
            )

        if self.backed_adata:
            indices = np.sort(indices)  # need to sort idxs for h5py datasets

        if np.amax(indices) >= len(self):
            raise IndexError("`indices` contains elements out of range.")

        sliced_data = {}
        for key, dtype in self.attributes_and_types.items():
            data = self.data[key]

            if isinstance(data, torch.Tensor):
                # if data is a torch tensor, we assume it is device backed
                data = data[indices]
            elif key == REGISTRY_KEYS.MINIFY_TYPE_KEY:
                # for minified  anndata, we need this because we can have a string
                # cur_data, which is the value of the MINIFY_TYPE_KEY in adata.uns,
                # used to record the type data minification
                # TODO: Adata manager should have a list of which fields it will load
                continue
            else:
                data = slice_and_convert(data, indices, dtype)

            # add singleton dimension since automatic batching is disabled
            if len(indices) == 1:
                data = data.reshape(1, -1)
            sliced_data[key] = data

        return sliced_data

    def __len__(self):
        return self.adata_manager.adata.shape[0]
