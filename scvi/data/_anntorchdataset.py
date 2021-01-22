import logging
from typing import Dict, List, Union

import anndata
import numpy as np
import pandas as pd
import torch
from torch.utils.data import Dataset

from scvi.data._anndata import get_from_registry

logger = logging.getLogger(__name__)


class AnnTorchDataset(Dataset):
    """Extension of torch dataset to get tensors from anndata."""

    def __init__(
        self,
        adata: anndata.AnnData,
        getitem_tensors: Union[List[str], Dict[str, type]] = None,
    ):
        self.adata = adata
        self.attributes_and_types = None
        self.getitem_tensors = getitem_tensors
        self.setup_getitem()
        self.setup_data_attr()

    @property
    def registered_keys(self):
        """Returns the keys of the mappings in scvi data registry."""
        return self.adata.uns["_scvi"]["data_registry"].keys()

    def setup_data_attr(self):
        """
        Sets data attribute.

        Reduces number of times anndata needs to be accessed
        """
        self.data = {
            key: get_from_registry(self.adata, key)
            for key, _ in self.attributes_and_types.items()
        }

    def setup_getitem(self):
        """
        Sets up the __getitem__ function used by Pytorch.

        By default, getitem will return every single item registered in the scvi data registry
        and will attempt to infer the correct type. np.float32 for continuous values, otherwise np.int64.

        If you want to specify which specific tensors to return you can pass in a List of keys from
        the scvi data registry. If you want to speficy specific tensors to return as well as their
        associated types, then you can pass in a dictionary with their type.

        Paramaters
        ----------
        getitem_tensors:
            Either a list of keys in the scvi data registry to return when getitem is called
            or

        Examples
        --------
        >>> sd = AnnTorchDataset(adata)

        # following will only return the X and batch_indices both by defualt as np.float32
        >>> sd.setup_getitem(getitem_tensors  = ['X,'batch_indices'])

        # This will return X as an integer and batch_indices as np.float32
        >>> sd.setup_getitem(getitem_tensors  = {'X':np.int64, 'batch_indices':np.float32])
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
            assert (
                key in registered_keys
            ), "{} not in anndata.uns['_scvi']['data_registry']".format(key)

        self.attributes_and_types = keys_to_type

    def __getitem__(self, idx: List[int]) -> Dict[str, torch.Tensor]:
        """Get tensors in dictionary from anndata at idx."""
        data_numpy = {}
        for key, dtype in self.attributes_and_types.items():
            data = self.data[key]
            if isinstance(data, np.ndarray):
                data_numpy[key] = data[idx].astype(dtype)
            elif isinstance(data, pd.DataFrame):
                data_numpy[key] = data.iloc[idx, :].to_numpy().astype(dtype)
            else:
                data_numpy[key] = data[idx].toarray().astype(dtype)

        return data_numpy

    def get_data(self, scvi_data_key):
        tensors = self.__getitem__(idx=[i for i in range(self.__len__())])
        return tensors[scvi_data_key]

    def __len__(self):
        return self.adata.shape[0]
