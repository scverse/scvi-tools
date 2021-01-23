import copy
import logging
from typing import Optional

import anndata
import numpy as np
import torch
from torch.utils.data import DataLoader

from scvi.data._anntorchdataset import AnnTorchDataset

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

    """

    def __init__(
        self,
        indices: np.ndarray,
        batch_size: int,
        shuffle: bool,
        drop_last: bool = False,
    ):
        self.indices = indices
        self.n_obs = len(indices)
        self.batch_size = batch_size
        self.shuffle = shuffle

    def __iter__(self):
        if self.shuffle is True:
            idx = torch.randperm(len(self.indices)).tolist()
        else:
            idx = torch.arange(len(self.indices)).tolist()

        data_iter = iter(
            [
                self.indices[idx[i : i + self.batch_size]]
                for i in range(0, len(idx), self.batch_size)
            ]
        )
        return data_iter

    def __len__(self):
        from math import ceil

        length = ceil(self.n_obs / self.batch_size)

        return length


class AnnDataLoader(DataLoader):
    """
    Scvi Data Loader for loading tensors from AnnData objects.

    Parameters
    ----------
    adata
        An anndata objects
    shuffle
        Whether the data should be shuffled
    indices
        The indices of the observations in the adata to load
    batch_size
        minibatch size to load each iteration
    data_and_attributes
        Dictionary with keys representing keys in data registry (`adata.uns["_scvi"]`)
        and value equal to desired numpy loading type (later made into torch tensor).
        If `None`, defaults to all registered data.
    data_loader_kwargs
        Keyword arguments for :class:`~torch.utils.data.DataLoader`
    """

    def __init__(
        self,
        adata: anndata.AnnData,
        shuffle=False,
        indices=None,
        batch_size=128,
        data_and_attributes: Optional[dict] = None,
        **data_loader_kwargs,
    ):

        if "_scvi" not in adata.uns.keys():
            raise ValueError("Please run setup_anndata() on your anndata object first.")

        if data_and_attributes is not None:
            data_registry = adata.uns["_scvi"]["data_registry"]
            for key in data_and_attributes.keys():
                if key not in data_registry.keys():
                    raise ValueError(
                        "{} required for model but not included when setup_anndata was run".format(
                            key
                        )
                    )

        self.dataset = AnnTorchDataset(adata, getitem_tensors=data_and_attributes)

        sampler_kwargs = {
            "batch_size": batch_size,
            "shuffle": shuffle,
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

        super().__init__(self.dataset, **self.data_loader_kwargs)
