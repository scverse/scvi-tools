import copy
import logging
from itertools import cycle
from typing import Dict, List, Optional, Sequence, Union

import numpy as np
import torch
from torch.utils.data import (
    BatchSampler,
    DataLoader,
    RandomSampler,
    SequentialSampler,
    Subset,
)

from scvi import REGISTRY_KEYS
from scvi.data import AnnDataManager
from scvi.data._utils import get_anndata_attribute

from ._datasets import AnnTorchDataset, DeviceBackedDataset

logger = logging.getLogger(__name__)


class AnnDataLoader(DataLoader):
    """DataLoader for loading tensors from AnnData objects.

    Parameters
    ----------
    adata_manager
        :class:`~scvi.data.AnnDataManager` object with a registered AnnData object.
    shuffle
        Whether the data should be shuffled.
    indices
        The indices of the observations in the adata to load.
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
        adata_manager: AnnDataManager,
        shuffle: bool = False,
        indices: Union[Sequence[int], Sequence[bool]] = None,
        batch_size: int = 128,
        data_and_attributes: Optional[Union[List[str], Dict[str, np.dtype]]] = None,
        drop_last: bool = False,
        iter_ndarray: bool = False,
        **data_loader_kwargs,
    ):
        self._full_dataset = AnnTorchDataset(
            adata_manager, getitem_tensors=data_and_attributes
        )
        if indices is None:
            indices = np.arange(len(self._full_dataset))
        else:
            if hasattr(indices, "dtype") and indices.dtype is np.dtype("bool"):
                indices = np.where(indices)[0].ravel()
            indices = np.asarray(indices)
        self.indices = indices

        # This is a lazy subset, it just remaps indices
        self.dataset = Subset(self._full_dataset, indices=self.indices)
        sampler_cls = SequentialSampler if not shuffle else RandomSampler
        sampler = BatchSampler(
            sampler=sampler_cls(self.dataset),
            batch_size=batch_size,
            drop_last=drop_last,
        )
        self.data_loader_kwargs = copy.deepcopy(data_loader_kwargs)
        # do not touch batch size here, sampler gives batched indices
        # This disables PyTorch automatic batching
        self.data_loader_kwargs.update({"sampler": sampler, "batch_size": None})

        if iter_ndarray:
            self.data_loader_kwargs.update({"collate_fn": _dummy_collate})

        super().__init__(self.dataset, **self.data_loader_kwargs)


def _dummy_collate(b: Dict[str, np.ndarray]) -> Dict[str, np.ndarray]:
    """Dummy collate to have dataloader return numpy ndarrays."""
    return b


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
        """Iter method for concat data loader.

        Will iter over the dataloader with the most data while cycling through
        the data in the other dataloaders. The order of data in returned iter_list
        is the same as indices_list.
        """
        iter_list = [
            cycle(dl) if dl != self.largest_dl else dl for dl in self.dataloaders
        ]
        return zip(*iter_list)


class SemiSupervisedDataLoader(ConcatDataLoader):
    """DataLoader that supports semisupervised training.

    Parameters
    ----------
    adata_manager
        :class:`~scvi.data.AnnDataManager` object that has been created via ``setup_anndata``.
    n_samples_per_label
        Number of subsamples for each label class to sample per epoch. By default, there
        is no label subsampling.
    indices
        The indices of the observations in the adata to load
    shuffle
        Whether the data should be shuffled
    batch_size
        minibatch size to load each iteration
    data_and_attributes
        Dictionary with keys representing keys in data registry (`adata_manager.data_registry`)
        and value equal to desired numpy loading type (later made into torch tensor).
        If `None`, defaults to all registered data.
    data_loader_kwargs
        Keyword arguments for :class:`~torch.utils.data.DataLoader`
    """

    def __init__(
        self,
        adata_manager: AnnDataManager,
        n_samples_per_label: Optional[int] = None,
        indices: Optional[List[int]] = None,
        shuffle: bool = False,
        batch_size: int = 128,
        data_and_attributes: Optional[dict] = None,
        drop_last: Union[bool, int] = False,
        **data_loader_kwargs,
    ):
        adata = adata_manager.adata
        if indices is None:
            indices = np.arange(adata.n_obs)

        self.indices = np.asarray(indices)

        if len(self.indices) == 0:
            return None

        self.n_samples_per_label = n_samples_per_label

        labels_state_registry = adata_manager.get_state_registry(
            REGISTRY_KEYS.LABELS_KEY
        )
        labels = get_anndata_attribute(
            adata_manager.adata,
            adata_manager.data_registry.labels.attr_name,
            labels_state_registry.original_key,
        ).ravel()

        # save a nested list of the indices per labeled category
        self.labeled_locs = []
        for label in np.unique(labels):
            if label != labels_state_registry.unlabeled_category:
                label_loc_idx = np.where(labels[indices] == label)[0]
                label_loc = self.indices[label_loc_idx]
                self.labeled_locs.append(label_loc)
        labelled_idx = self.subsample_labels()

        super().__init__(
            adata_manager=adata_manager,
            indices_list=[self.indices, labelled_idx],
            shuffle=shuffle,
            batch_size=batch_size,
            data_and_attributes=data_and_attributes,
            drop_last=drop_last,
            **data_loader_kwargs,
        )

    def resample_labels(self):
        """Resamples the labeled data."""
        labelled_idx = self.subsample_labels()
        # self.dataloaders[0] iterates over full_indices
        # self.dataloaders[1] iterates over the labelled_indices
        # change the indices of the labelled set
        self.dataloaders[1].indices = labelled_idx

    def subsample_labels(self):
        """Subsamples each label class by taking up to n_samples_per_label samples per class."""
        if self.n_samples_per_label is None:
            return np.concatenate(self.labeled_locs)

        sample_idx = []
        for loc in self.labeled_locs:
            if len(loc) < self.n_samples_per_label:
                sample_idx.append(loc)
            else:
                label_subset = np.random.choice(
                    loc, self.n_samples_per_label, replace=False
                )
                sample_idx.append(label_subset)
        sample_idx = np.concatenate(sample_idx)
        return sample_idx

    class DeviceBackedDataLoader(DataLoader):
        """DataLoader for loading device-backed tensors."""

        def __init__(
            self,
            adata_manager: AnnDataManager,
            device: torch.device,
            shuffle: bool = False,
            indices: Union[Sequence[int], Sequence[bool]] = None,
            batch_size: int = 128,
            data_and_attributes: Optional[Union[List[str], Dict[str, np.dtype]]] = None,
            drop_last: bool = False,
            iter_ndarray: bool = False,
            **kwargs,
        ):
            tensor_dict = self._get_tensor_dict(
                adata_manager, indices, device, **kwargs
            )
            dataset = DeviceBackedDataset(tensor_dict)
            batch_size = batch_size or len(dataset)
            sampler_cls = RandomSampler if shuffle else SequentialSampler
            sampler = BatchSampler(
                sampler=sampler_cls(dataset),
                batch_size=batch_size,
                drop_last=kwargs.pop("drop_last", False),
            )
            super().__init__(dataset, sampler=sampler, batch_size=None)

        def _get_tensor_dict(
            adata_manager: AnnDataManager,
            indices: Union[Sequence[int], Sequence[bool]],
            device: torch.device,
            **kwargs,
        ) -> Dict[str, torch.Tensor]:
            """Get tensor dict for a given set of indices."""
            if len(indices) is None or len(indices) == 0:
                return

            dl = AnnDataLoader(
                adata_manager,
                indices=indices,
                batch_size=len(indices),
                shuffle=False,
                pin_memory=kwargs.pop("pin_memory", False),
                **kwargs,
            )
            batch = next(iter(dl))

            return {k: v.to(device) for k, v in batch.items()}
