import logging
from typing import Dict, List, Union

import anndata
import numpy as np
import pandas as pd
import torch
from torch.utils.data import Dataset

from scvi.data._anndata import get_from_registry

logger = logging.getLogger(__name__)


class ScviDataset(Dataset):
    """Extension of torch dataset to get tensors from anndata."""

    def __init__(
        self,
        adata: anndata.AnnData,
        getitem_tensors: Union[List[str], Dict[str, type]] = None,
    ):
        self.adata = adata
        self.attributes_and_types = None
        self.setup_getitem(getitem_tensors)
        self.setup_data_attr()
        self.gene_names = self.adata.var_names
        self.normalized_X = None

    def get_registered_keys(
        self,
    ):
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

    def setup_getitem(self, getitem_tensors: Union[List[str], Dict[str, type]] = None):
        """
        Sets up the __getitem__ function used by Pytorch.

        By default, getitem will return every single item registered in the scvi data registry
        and also set their type to np.float32.

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
        >>> sd = ScviDataset(adata)

        # following will only return the X and batch_indices both by defualt as np.float32
        >>> sd .setup_getitem(getitem_tensors  = ['X,'batch_indices'])

        # This will return X as an integer and batch_indices as np.float32
        >>> sd.setup_getitem(getitem_tensors  = {'X':np.int64, 'batch_indices':np.float32])
        """
        registered_keys = self.get_registered_keys()

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
            ), "{} not in anndata.uns['scvi_data_registry']".format(key)

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

    @property
    def n_cells(self) -> int:
        """Returns the number of cells in the dataset."""
        n_cells = self.adata.uns["_scvi"]["summary_stats"]["n_cells"]
        return n_cells

    @property
    def n_vars(self) -> int:
        """Returns the number of variables in the dataset."""
        n_vars = self.adata.uns["_scvi"]["summary_stats"]["n_vars"]
        return n_vars

    @property
    def n_batches(self) -> int:
        return self.adata.uns["_scvi"]["summary_stats"]["n_batch"]

    def to_anndata(
        self,
    ) -> anndata.AnnData:
        return self.adata.copy()
