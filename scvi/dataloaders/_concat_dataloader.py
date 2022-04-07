from itertools import cycle
from typing import List, Optional, Union

import numpy as np
from torch.utils.data import DataLoader

from scvi.data import AnnDataManager

from ._ann_dataloader import AnnDataLoader


class ConcatDataLoader(DataLoader):
    """
    DataLoader that supports a list of list of indices to load.

    Parameters
    ----------
    adata_manager
        :class:`~scvi.data.AnnDataManager` object that has been created via ``setup_anndata``.
    indices_list
        List where each element is a list of indices in the adata to load
    shuffle
        Whether the data should be shuffled
    batch_size
        minibatch size to load each iteration
    data_and_attributes
        Dictionary with keys representing keys in data registry (``adata_manager.data_registry``)
        and value equal to desired numpy loading type (later made into torch tensor).
        If ``None``, defaults to all registered data.
    data_loader_kwargs
        Keyword arguments for :class:`~torch.utils.data.DataLoader`
    """

    def __init__(
        self,
        adata_manager: AnnDataManager,
        indices_list: List[List[int]],
        shuffle: bool = False,
        batch_size: int = 128,
        data_and_attributes: Optional[dict] = None,
        drop_last: Union[bool, int] = False,
        **data_loader_kwargs,
    ):
        self.dataloaders = []
        for indices in indices_list:
            self.dataloaders.append(
                AnnDataLoader(
                    adata_manager,
                    indices=indices,
                    shuffle=shuffle,
                    batch_size=batch_size,
                    data_and_attributes=data_and_attributes,
                    drop_last=drop_last,
                    **data_loader_kwargs,
                )
            )
        lens = [len(dl) for dl in self.dataloaders]
        self.largest_dl = self.dataloaders[np.argmax(lens)]
        super().__init__(self.largest_dl, **data_loader_kwargs)

    def __len__(self):
        return len(self.largest_dl)

    def __iter__(self):
        """
        Iter method for concat data loader.

        Will iter over the dataloader with the most data while cycling through
        the data in the other dataloaders. The order of data in returned iter_list
        is the same as indices_list.
        """
        iter_list = [
            cycle(dl) if dl != self.largest_dl else dl for dl in self.dataloaders
        ]
        return zip(*iter_list)
