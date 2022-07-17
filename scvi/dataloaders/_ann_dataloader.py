import copy
import logging
from typing import Optional, Union

import numpy as np
from catalyst.data.sampler import DistributedSamplerWrapper
from torch.utils.data import DataLoader

from scvi.data import AnnDataManager

from ._anntorchdataset import AnnTorchDataset
from ._samplers import BatchSampler, SubsetDistributedSampler

logger = logging.getLogger(__name__)


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
        sampler_cls: Union[BatchSampler, SubsetDistributedSampler] = BatchSampler,
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
            "shuffle": shuffle,
            "drop_last": drop_last,
        }
        sampler_kwargs.update({"batch_size": batch_size})

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

        try:
            import torch_xla.core.xla_model as xm

            distributed_sampler_kwargs = dict(
                num_replicas=xm.xrt_world_size(), rank=xm.get_ordinal()
            )
        except ModuleNotFoundError:
            distributed_sampler_kwargs = {}

        sampler = BatchSampler(**self.sampler_kwargs)
        self.data_loader_kwargs = copy.copy(data_loader_kwargs)

        if sampler_cls == SubsetDistributedSampler:
            sampler = DistributedSamplerWrapper(sampler, **distributed_sampler_kwargs)
        # do not touch batch size here, sampler gives batched indices
        self.data_loader_kwargs.update({"sampler": sampler, "batch_size": None})

        if iter_ndarray:
            self.data_loader_kwargs.update({"collate_fn": _dummy_collate})

        super().__init__(self.dataset, **self.data_loader_kwargs)


def _dummy_collate(b):
    """Dummy collate to have dataloader return numpy ndarrays."""
    return b
