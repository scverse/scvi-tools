import copy
import logging
from typing import Dict, List, Optional, Sequence, Union

import numpy as np
from torch.utils.data import (
    BatchSampler,
    DataLoader,
    Dataset,
    RandomSampler,
    Sampler,
    SequentialSampler,
    Subset,
)

from scvi.data import AnnDataManager

from ._anntorchdataset import AnnTorchDataset

logger = logging.getLogger(__name__)


class AnnDataLoader(DataLoader):
    """DataLoader for loading tensors from AnnData objects.

    Parameters
    ----------
    adata_manager
        :class:`~scvi.data.AnnDataManager` object that has been created via a model's
        `setup_anndata` or `setup_mudata` method. Must be provided if `dataset` is not.
    dataset
        The :class:`~torch.utils.data.Dataset` to load from. Must be provided if
        `adata_manager` is not or if `sampler` is provided.
    shuffle
        Whether the data should be shuffled.
    indices
        The indices of the observations in the adata to load.
    sampler
        The sampler used for the dataloader. If `None`, a
        :class:`~torch.utils.data.BatchSampler` will be used.
    batch_size
        minibatch size to load each iteration.
    data_and_attributes
        Dictionary with keys representing keys in data registry (``adata_manager.data_registry``)
        and value equal to desired numpy loading type (later made into torch tensor) or list of
        such keys. A list can be used to subset to certain keys in the event that more tensors than
        needed have been registered. If ``None``, defaults to all registered data.
    iter_ndarray
        Whether to iterate over numpy arrays instead of torch tensors.
    data_loader_kwargs
        Keyword arguments for :class:`~torch.utils.data.DataLoader`.
    """

    def __init__(
        self,
        adata_manager: Optional[AnnDataManager] = None,
        dataset: Optional[Dataset] = None,
        shuffle: bool = False,
        indices: Union[Sequence[int], Sequence[bool]] = None,
        sampler: Optional[Sampler] = None,
        batch_size: int = 128,
        data_and_attributes: Optional[Union[List[str], Dict[str, np.dtype]]] = None,
        drop_last: bool = False,
        iter_ndarray: bool = False,
        **data_loader_kwargs,
    ):
        if adata_manager is not None and dataset is not None:
            raise ValueError("Cannot provide both `adata_manager` and `dataset`.")
        if adata_manager is None and dataset is None:
            raise ValueError("Must provide either `adata_manager` or `dataset`.")
        if sampler is not None and dataset is None:
            raise ValueError("Must provide `dataset` if `sampler` is provided.")

        if dataset is None:
            self._full_dataset = AnnTorchDataset(
                adata_manager, getitem_tensors=data_and_attributes
            )
        else:
            self._full_dataset = dataset

        if indices is None:
            indices = np.arange(len(self._full_dataset))
        else:
            if hasattr(indices, "dtype") and indices.dtype is np.dtype("bool"):
                indices = np.where(indices)[0].ravel()
            indices = np.asarray(indices)
        self.indices = indices

        # This is a lazy subset, it just remaps indices
        self.dataset = Subset(self._full_dataset, indices=self.indices)

        if sampler is None:
            sampler_cls = SequentialSampler if not shuffle else RandomSampler
            sampler = BatchSampler(
                sampler=sampler_cls(self.dataset),
                batch_size=batch_size,
                drop_last=drop_last,
            )
        if isinstance(sampler, BatchSampler):
            batch_size = None  # disables PyTorch automatic batching
            shuffle = False
            drop_last = False

        self.data_loader_kwargs = copy.deepcopy(data_loader_kwargs)

        if iter_ndarray:
            self.data_loader_kwargs.update({"collate_fn": _dummy_collate})

        super().__init__(
            self.dataset,
            sampler=sampler,
            batch_size=batch_size,
            shuffle=shuffle,
            drop_last=drop_last,
            **self.data_loader_kwargs,
        )


def _dummy_collate(b: Dict[str, np.ndarray]) -> Dict[str, np.ndarray]:
    """Dummy collate to have dataloader return numpy ndarrays."""
    return b
