from itertools import cycle

import numpy as np
from torch.utils.data import DataLoader

from scvi.data import AnnDataManager

from ._ann_dataloader import AnnDataLoader


class ConcatDataLoader(DataLoader):
    """DataLoader that supports a list of list of indices to load.

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
    drop_last
        If `True` and the dataset is not evenly divisible by `batch_size`, the last
        incomplete batch is dropped. If `False` and the dataset is not evenly divisible
        by `batch_size`, then the last batch will be smaller than `batch_size`.
    distributed_sampler
        ``EXPERIMENTAL`` Whether to use :class:`~scvi.dataloaders.BatchDistributedSampler` as the
        sampler. If `True`, `sampler` must be `None`.
    data_loader_kwargs
        Keyword arguments for :class:`~torch.utils.data.DataLoader`
    """

    def __init__(
        self,
        adata_manager: AnnDataManager,
        indices_list: list[list[int]],
        shuffle: bool = False,
        batch_size: int = 128,
        data_and_attributes: dict | None = None,
        drop_last: bool | int = False,
        distributed_sampler: bool = False,
        **data_loader_kwargs,
    ):
        self.adata_manager = adata_manager
        self.data_loader_kwargs = data_loader_kwargs
        self.data_and_attributes = data_and_attributes
        self._shuffle = shuffle
        self._batch_size = batch_size
        self._drop_last = drop_last
        self._distributed_sampler = distributed_sampler

        self.dataloaders = []
        for indices in indices_list:
            if self._distributed_sampler:
                self.data_loader_kwargs.pop("sampler", None)
            self.dataloaders.append(
                AnnDataLoader(
                    adata_manager,
                    indices=indices,
                    shuffle=shuffle,
                    batch_size=batch_size,
                    data_and_attributes=data_and_attributes,
                    drop_last=drop_last,
                    distributed_sampler=distributed_sampler,
                    **self.data_loader_kwargs,
                )
            )
        lens = [len(dl) for dl in self.dataloaders]
        self.largest_dl = self.dataloaders[np.argmax(lens)]
        self.data_loader_kwargs.pop("drop_dataset_tail", None)
        super().__init__(self.largest_dl, **self.data_loader_kwargs)

    def __len__(self):
        return len(self.largest_dl)

    def __iter__(self):
        """Iter method for concat data loader.

        Will iter over the dataloader with the most data while cycling through
        the data in the other dataloaders. The order of data in returned iter_list
        is the same as indices_list.
        """
        if not self._distributed_sampler:
            iter_list = [cycle(dl) if dl != self.largest_dl else dl for dl in self.dataloaders]
            strict = True
        else:
            iter_list = self.dataloaders
            strict = False
        return zip(*iter_list, strict=strict)
