import logging
import warnings
from itertools import cycle
from typing import Dict, List, Optional, Sequence, Union

import numpy as np
from torch.utils.data import (
    BatchSampler,
    DataLoader,
    RandomSampler,
    Sampler,
    SequentialSampler,
    Subset,
)

from scvi import REGISTRY_KEYS, settings
from scvi.data import AnnDataManager
from scvi.data._utils import get_anndata_attribute

from ._datasets import AnnTorchDataset
from ._docstrings import dataloader_dsp

logger = logging.getLogger(__name__)


@dataloader_dsp.dedent
class AnnDataLoader(DataLoader):
    """%(summary)s

    Parameters
    ----------
    %(param_adata_manager)s
    %(param_indices)s
    %(param_sampler)s
    %(param_shuffle)s
    %(param_batch_size)s
    %(param_data_and_attributes)s
    %(param_drop_last)s
    %(param_iter_ndarray)s
    %(param_pin_memory)s
    %(param_accelerator)s
    %(param_device)s
    %(param_device_backed)s
    %(param_kwargs)s
    """

    def __init__(
        self,
        adata_manager: AnnDataManager,
        indices: Union[Sequence[int], Sequence[bool]] = None,
        sampler: Optional[Sampler] = None,
        shuffle: bool = False,
        batch_size: int = 128,
        data_and_attributes: Optional[Union[List[str], Dict[str, np.dtype]]] = None,
        drop_last: bool = False,
        iter_ndarray: bool = False,
        pin_memory: bool = False,
        accelerator: str = "auto",
        device: Union[int, str] = "auto",
        device_backed: bool = False,
        **kwargs,
    ):
        dataset = AnnTorchDataset(
            adata_manager,
            getitem_tensors=data_and_attributes,
            accelerator=accelerator,
            device=device,
            device_backed=device_backed,
        )
        if indices is None:
            indices = np.arange(len(dataset))
        else:
            if hasattr(indices, "dtype") and indices.dtype is np.dtype("bool"):
                indices = np.where(indices)[0].ravel()
            indices = np.asarray(indices)

        if pin_memory and device_backed:
            warnings.warn(
                "Cannot set `pin_memory=True` when `device_backed=True`. "
                "Setting `pin_memory=False`.",
                UserWarning,
                stacklevel=settings.warnings_stacklevel,
            )

        if iter_ndarray:
            if "collate_fn" in kwargs:
                raise ValueError("Cannot set `collate_fn` when `iter_ndarray=True`. ")
            kwargs["collate_fn"] = lambda x: x

        dataset = Subset(dataset, indices=indices)  # lazy subset, remaps indices

        if sampler is None:
            sampler = RandomSampler(dataset) if shuffle else SequentialSampler(dataset)
            sampler = BatchSampler(
                sampler=sampler,
                batch_size=batch_size,
                drop_last=drop_last,
            )
            batch_size = None  # disables torch automatic batching
            shuffle = False
            drop_last = False

        super().__init__(
            dataset,
            sampler=sampler,
            batch_size=batch_size,
            shuffle=shuffle,
            drop_last=drop_last,
            pin_memory=pin_memory,
            **kwargs,
        )


@dataloader_dsp.dedent
class ConcatDataLoader(DataLoader):
    """%(summary)s

    Supports loading a list of list of indices.

    Parameters
    ----------
    %(param_adata_manager)s
    %(param_indices_list)s
    %(param_shuffle)s
    %(param_batch_size)s
    %(param_data_and_attributes)s
    %(param_drop_last)s
    %(param_iter_ndarray)s
    %(param_accelerator)s
    %(param_device)s
    %(param_device_backed)s
    %(param_kwargs)s
    """

    def __init__(
        self,
        adata_manager: AnnDataManager,
        indices_list: List[List[int]],
        shuffle: bool = False,
        batch_size: int = 128,
        data_and_attributes: Optional[dict] = None,
        drop_last: Union[bool, int] = False,
        iter_ndarray: bool = False,
        accelerator: str = "auto",
        device: Union[int, str] = "auto",
        device_backed: bool = False,
        **kwargs,
    ):
        self.dataloaders = []
        n_obs_dls = []

        for indices in indices_list:
            dl = AnnDataLoader(
                adata_manager,
                indices=indices,
                shuffle=shuffle,
                batch_size=batch_size,
                data_and_attributes=data_and_attributes,
                drop_last=drop_last,
                iter_ndarray=iter_ndarray,
                accelerator=accelerator,
                device=device,
                device_backed=device_backed,
                **kwargs,
            )
            self.dataloaders.append(dl)
            n_obs_dls.append(len(dl))

        self.largest_dl = self.dataloaders[np.argmax(n_obs_dls)]
        super().__init__(self.largest_dl, **kwargs)

    def __len__(self):
        return len(self.largest_dl)

    def __iter__(self):
        """Iterate over the dataloaders.

        Iterates over the dataloader with the most data while cycling through the data
        in the other datalaoders. The order of data returned is the same as that
        passed in `indices_list`.
        """
        iter_list = [
            cycle(dl) if dl != self.largest_dl else dl for dl in self.dataloaders
        ]
        return zip(*iter_list)


@dataloader_dsp.dedent
class SemiSupervisedDataLoader(ConcatDataLoader):
    """%(summary)s

    Supports loading labeled and unlabeled data.

    Parameters
    ----------
    %(param_adata_manager)s
    %(param_indices)s
    %(param_n_samples_per_label)s
    %(param_shuffle)s
    %(param_batch_size)s
    %(param_data_and_attributes)s
    %(param_drop_last)s
    %(param_iter_ndarray)s
    %(param_accelerator)s
    %(param_device)s
    %(param_device_backed)s
    %(param_kwargs)s
    """

    def __init__(
        self,
        adata_manager: AnnDataManager,
        indices: Optional[List[int]] = None,
        n_samples_per_label: Optional[int] = None,
        shuffle: bool = False,
        batch_size: int = 128,
        data_and_attributes: Optional[dict] = None,
        drop_last: Union[bool, int] = False,
        iter_ndarray: bool = False,
        accelerator: str = "auto",
        device: Union[int, str] = "auto",
        device_backed: bool = False,
        **kwargs,
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
            iter_ndarray=iter_ndarray,
            accelerator=accelerator,
            device=device,
            device_backed=device_backed,
            **kwargs,
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
