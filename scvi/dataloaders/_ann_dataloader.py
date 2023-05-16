import copy
import logging
from typing import Dict, List, Optional, Sequence, Union

import numpy as np
from torch.utils.data import (
    BatchSampler,
    DataLoader,
    RandomSampler,
    Sampler,
    SequentialSampler,
)

from scvi.data import AnnDataManager

logger = logging.getLogger(__name__)


class AnnDataLoader(DataLoader):
    """DataLoader for loading tensors from AnnData objects.

    Parameters
    ----------
    adata_manager
        :class:`~scvi.data.AnnDataManager` object with a registered AnnData object.
    shuffle
        Whether the data should be shuffled
    indices
        The indices of the observations in the adata to load
    batch_size
        minibatch size to load each iteration
    sampler
        Defines the strategy to draw samples from the dataset. Can be any Iterable with __len__ implemented.
        If specified, shuffle must not be specified. By default, we use a custom sampler that is designed to
        get a minibatch of data with one call to __getitem__.
    data_and_attributes
        Dictionary with keys representing keys in data registry (``adata_manager.data_registry``)
        and value equal to desired numpy loading type (later made into torch tensor) or list of
        such keys. A list can be used to subset to certain keys in the event that more tensors than
        needed have been registered. If ``None``, defaults to all registered data.
    iter_ndarray
        Whether to iterate over numpy arrays instead of torch tensors
    data_loader_kwargs
        Keyword arguments for :class:`~torch.utils.data.DataLoader`
    """

    def __init__(
        self,
        adata_manager: AnnDataManager,
        shuffle: bool = False,
        indices: Union[Sequence[int], Sequence[bool]] = None,
        batch_size: int = 128,
        sampler: Optional[Sampler] = None,
        data_and_attributes: Optional[Union[List[str], Dict[str, np.dtype]]] = None,
        drop_last: bool = False,
        iter_ndarray: bool = False,
        **data_loader_kwargs,
    ):
        if indices is None:
            indices = np.arange(adata_manager.adata.shape[0])
        else:
            if hasattr(indices, "dtype") and indices.dtype is np.dtype("bool"):
                indices = np.where(indices)[0].ravel()
            indices = np.asarray(indices)
        self.indices = indices
        self.dataset = adata_manager.create_torch_dataset(
            indices=indices, data_and_attributes=data_and_attributes
        )
        self.data_loader_kwargs = copy.deepcopy(data_loader_kwargs)
        if sampler is None:
            sampler_cls = SequentialSampler if not shuffle else RandomSampler
            sampler = BatchSampler(
                sampler=sampler_cls(self.dataset),
                batch_size=batch_size,
                drop_last=drop_last,
            )
            # do not touch batch size here, sampler gives batched indices
            # This disables PyTorch automatic batching, which is necessary
            # for fast access to sparse matrices
            self.data_loader_kwargs.update({"sampler": sampler, "batch_size": None})

        if iter_ndarray:
            self.data_loader_kwargs.update({"collate_fn": _dummy_collate})

        super().__init__(self.dataset, **self.data_loader_kwargs)


def _dummy_collate(b: Dict[str, np.ndarray]) -> Dict[str, np.ndarray]:
    """Dummy collate to have dataloader return numpy ndarrays."""
    return b
