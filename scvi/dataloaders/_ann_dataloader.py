import copy
import logging
from typing import Optional, Union

import numpy as np
import torch
from torch.utils.data import DataLoader

from scvi.data import AnnDataManager

from ._anntorchdataset import AnnTorchDataset

logger = logging.getLogger(__name__)


class BatchSampler(torch.utils.data.sampler.Sampler):
    """
    Custom torch Sampler that returns a list of indices of size batch_size.

    Parameters
    ----------
    indices
        list of indices to sample from
    batch_size
        batch size of each iteration
    shuffle
        if ``True``, shuffles indices before sampling
    drop_last
        if int, drops the last batch if its length is less than drop_last.
        if drop_last == True, drops last non-full batch.
        if drop_last == False, iterate over all batches.
    """

    def __init__(
        self,
        indices: np.ndarray,
        batch_size: int,
        shuffle: bool,
        drop_last: Union[bool, int] = False,
    ):
        self.indices = indices
        self.n_obs = len(indices)
        self.batch_size = batch_size
        self.shuffle = shuffle

        if drop_last > batch_size:
            raise ValueError(
                "drop_last can't be greater than batch_size. "
                + f"drop_last is {drop_last} but batch_size is {batch_size}."
            )

        last_batch_len = self.n_obs % self.batch_size
        if (drop_last is True) or (last_batch_len < drop_last):
            drop_last_n = last_batch_len
        elif (drop_last is False) or (last_batch_len >= drop_last):
            drop_last_n = 0
        else:
            raise ValueError("Invalid input for drop_last param. Must be bool or int.")

        self.drop_last_n = drop_last_n

    def __iter__(self):
        if self.shuffle is True:
            idx = torch.randperm(self.n_obs).tolist()
        else:
            idx = torch.arange(self.n_obs).tolist()

        if self.drop_last_n != 0:
            idx = idx[: -self.drop_last_n]

        data_iter = iter(
            [
                self.indices[idx[i : i + self.batch_size]]
                for i in range(0, len(idx), self.batch_size)
            ]
        )
        return data_iter

    def __len__(self):
        from math import ceil

        if self.drop_last_n != 0:
            length = self.n_obs // self.batch_size
        else:
            length = ceil(self.n_obs / self.batch_size)
        return length


class AnnDataLoader(DataLoader):
    """
    DataLoader for loading tensors from AnnData objects.

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
    data_and_attributes
        Dictionary with keys representing keys in data registry (``adata_manager.data_registry``)
        and value equal to desired numpy loading type (later made into torch tensor).
        If ``None``, defaults to all registered data.
    data_loader_kwargs
        Keyword arguments for :class:`~torch.utils.data.DataLoader`
    iter_ndarray
        Whether to iterate over numpy arrays instead of torch tensors
    """

    def __init__(
        self,
        adata_manager: AnnDataManager,
        shuffle=False,
        indices=None,
        batch_size=128,
        data_and_attributes: Optional[dict] = None,
        drop_last: Union[bool, int] = False,
        iter_ndarray: bool = False,
        **data_loader_kwargs,
    ):

        if adata_manager.adata is None:
            raise ValueError(
                "Please run register_fields() on your AnnDataManager object first."
            )

        if data_and_attributes is not None:
            data_registry = adata_manager.data_registry
            for key in data_and_attributes.keys():
                if key not in data_registry.keys():
                    raise ValueError(
                        f"{key} required for model but not registered with AnnDataManager."
                    )

        self.dataset = AnnTorchDataset(
            adata_manager, getitem_tensors=data_and_attributes
        )

        sampler_kwargs = {
            "batch_size": batch_size,
            "shuffle": shuffle,
            "drop_last": drop_last,
        }

        if indices is None:
            indices = np.arange(len(self.dataset))
            sampler_kwargs["indices"] = indices
        else:
            if hasattr(indices, "dtype") and indices.dtype is np.dtype("bool"):
                indices = np.where(indices)[0].ravel()
            indices = np.asarray(indices)
            sampler_kwargs["indices"] = indices

        self.indices = indices
        self.sampler_kwargs = sampler_kwargs
        sampler = BatchSampler(**self.sampler_kwargs)
        self.data_loader_kwargs = copy.copy(data_loader_kwargs)
        # do not touch batch size here, sampler gives batched indices
        self.data_loader_kwargs.update({"sampler": sampler, "batch_size": None})

        if iter_ndarray:
            self.data_loader_kwargs.update({"collate_fn": _dummy_collate})

        super().__init__(self.dataset, **self.data_loader_kwargs)


def _dummy_collate(b):
    """Dummy collate to have dataloader return numpy ndarrays."""
    return b
