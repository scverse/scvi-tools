import logging
import warnings
from itertools import cycle
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
    %(param_dataset)s
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
        adata_manager: Optional[AnnDataManager] = None,
        dataset: Optional[Dataset] = None,
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
        if adata_manager is not None and dataset is not None:
            raise ValueError("Cannot set both `adata_manager` and `dataset`.")
        if adata_manager is None and dataset is None:
            raise ValueError("Must set either `adata_manager` or `dataset`.")
        if adata_manager is not None:
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
            if len(indices) == 0:
                raise ValueError("`indices` must contain at least one index.")

        if pin_memory and device_backed:
            warnings.warn(
                "Cannot set `pin_memory=True` when `device_backed=True`. "
                "Setting `pin_memory=False`.",
                UserWarning,
                stacklevel=settings.warnings_stacklevel,
            )
            pin_memory = False

        if iter_ndarray:
            if device_backed:
                raise ValueError(
                    "Cannot set `iter_ndarray=True` when `device_backed=True`."
                )
            if "collate_fn" in kwargs:
                raise ValueError("Cannot set `collate_fn` when `iter_ndarray=True`.")
            kwargs["collate_fn"] = lambda x: x

        dataset = Subset(dataset, indices=indices)  # lazy subset, remaps indices

        if sampler is None:
            sampler = RandomSampler(dataset) if shuffle else SequentialSampler(dataset)
            sampler = BatchSampler(
                sampler=sampler,
                batch_size=batch_size,
                drop_last=drop_last,
            )
            # disables torch automatic batching
            batch_size = None
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

    Supports loading a list of list of indices. Iterates over the largest indices list
    while cycling through the data in the other lists.

    Parameters
    ----------
    %(param_adata_manager)s
    %(param_indices_list)s
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
        indices_list: List[List[int]],
        sampler: Optional[Sampler] = None,
        shuffle: bool = False,
        batch_size: int = 128,
        data_and_attributes: Optional[dict] = None,
        drop_last: Union[bool, int] = False,
        iter_ndarray: bool = False,
        pin_memory: bool = False,
        accelerator: str = "auto",
        device: Union[int, str] = "auto",
        device_backed: bool = False,
        **kwargs,
    ):
        self._adata_manager = adata_manager
        self._dataset = AnnTorchDataset(
            adata_manager,
            getitem_tensors=data_and_attributes,
            accelerator=accelerator,
            device=device,
            device_backed=device_backed,
        )
        self._dataloader_kwargs = {
            "sampler": sampler,
            "shuffle": shuffle,
            "batch_size": batch_size,
            "data_and_attributes": data_and_attributes,
            "drop_last": drop_last,
            "iter_ndarray": iter_ndarray,
            "pin_memory": pin_memory,
            "accelerator": accelerator,
            "device": device,
            "device_backed": device_backed,
        }
        self._dataloader_kwargs.update(kwargs)
        self.indices_list = indices_list

    @property
    def dataloaders(self) -> List[DataLoader]:
        """Returns the list of dataloaders."""
        return self._dataloaders

    @dataloaders.setter
    def dataloaders(self, value: List[DataLoader]):
        """Sets the list of dataloaders."""
        raise AttributeError(
            "Cannot set `dataloaders` after initialization. Please set "
            "`indices_list` instead."
        )

    @property
    def indices_list(self) -> list:
        """Returns the list of indices."""
        return self._indices_list

    @indices_list.setter
    def indices_list(self, value: list):
        """Sets the list of indices."""
        dataloaders = []
        for indices in value:
            # passing in the same dataset to all dataloaders ensures that if there are
            # non-exclusive indices, observations are only loaded once and shared
            # between dataloaders
            dataloader = AnnDataLoader(
                dataset=self._dataset, indices=indices, **self._dataloader_kwargs
            )
            dataloaders.append(dataloader)

        self._indices_list = value
        self._dataloaders = dataloaders
        _ = len(self)  # recompute length

    def __len__(self):
        """Returns the length of the largest dataloader."""
        if hasattr(self, "_length"):
            return self._length

        self._length = max([len(dl) for dl in self.dataloaders])
        return self._length

    def __iter__(self):
        """Iterate over the dataloaders.

        Iterates over the dataloader with the most data while cycling through the data
        in the other datalaoders. The order of data returned is the same as that
        passed in `indices_list`.
        """
        iter_order = [
            cycle(dl) if len(dl) < len(self) else dl for dl in self.dataloaders
        ]
        return zip(*iter_order)


@dataloader_dsp.dedent
class SemiSupervisedDataLoader(ConcatDataLoader):
    """%(summary)s

    Labels must be present and registered with :class:`scvi.data.AnnDataManager`
    to use this dataloader. Loads two minibatches, one that contains both
    labeled and unlabeled data, and one that contains only labeled data.

    Parameters
    ----------
    %(param_adata_manager)s
    %(param_indices)s
    %(param_n_samples_per_label)s
    %(param_sampler)s
    %(param_shuffle)s
    %(param_batch_size)s
    %(param_data_and_attributes)s
    %(param_drop_last)s
    %(param_iter_ndarray)s
    %(param_seed)s
    %(param_pin_memory)s
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
        sampler: Optional[Sampler] = None,
        shuffle: bool = False,
        batch_size: int = 128,
        data_and_attributes: Optional[dict] = None,
        drop_last: Union[bool, int] = False,
        iter_ndarray: bool = False,
        seed: Optional[int] = None,
        pin_memory: bool = False,
        accelerator: str = "auto",
        device: Union[int, str] = "auto",
        device_backed: bool = False,
        **kwargs,
    ):
        if REGISTRY_KEYS.LABELS_KEY not in adata_manager.data_registry:
            raise ValueError(
                "Must have labels in the data registry for `SemiSupervisedDataLoader`."
            )

        if indices is None:
            indices = np.arange(adata_manager.adata.n_obs)
        seed = seed or settings.seed

        super().__init__(
            adata_manager,
            indices_list=[indices],  # dummy indices_list as it needs to be set later
            sampler=sampler,
            shuffle=shuffle,
            batch_size=batch_size,
            data_and_attributes=data_and_attributes,
            drop_last=drop_last,
            iter_ndarray=iter_ndarray,
            pin_memory=pin_memory,
            accelerator=accelerator,
            device=device,
            device_backed=device_backed,
            **kwargs,
        )

        self.indices = indices
        self.generator = np.random.default_rng(seed=seed)
        self.n_samples_per_label = n_samples_per_label
        self.resample_labels()

    @property
    def n_samples_per_label(self) -> int:
        """Returns the number of samples per label."""
        return self._n_samples_per_label

    @n_samples_per_label.setter
    def n_samples_per_label(self, value: int):
        if hasattr(self, "_n_samples_per_label") and self._n_samples_per_label != value:
            # if the value is changed, resample the labels
            self.resample_labels()

        if value is not None:
            min_n_samples = min(self.labeled_indices.values())
            if value > min_n_samples:
                warnings.warn(
                    "n_samples_per_label` is greater than the number of samples in the "
                    "smallest label. Labels with less than `n_samples_per_label` will "
                    "be sampled with replacement.",
                    UserWarning,
                    stacklevel=settings.warnings_stacklevel,
                )
        self._n_samples_per_label = value

    @property
    def labeled_indices(self) -> dict:
        if hasattr(self, "_labeled_indices"):
            return self._labeled_indices

        labels_state_registry = self._adata_manager.get_state_registry(
            REGISTRY_KEYS.LABELS_KEY
        )
        labels = get_anndata_attribute(
            self._adata_manager.adata,
            self._adata_manager.data_registry.labels.attr_name,
            labels_state_registry.original_key,
        ).ravel()[self.indices]

        indices = {}
        for label in np.unique(labels):
            if hasattr(labels_state_registry, "unlabeled_category") and (
                label == labels_state_registry.unlabeled_category
            ):
                continue
            locs = self.indices[np.where(labels == label)[0]]
            indices[label] = locs

        self._labeled_indices = indices
        return self._labeled_indices

    @property
    def subsampled_labeled_indices(self) -> List[int]:
        """Resamples and returns the subsampled labeled indices."""
        if hasattr(self, "_subsampled_labeled_indices"):
            return self._subsampled_labeled_indices

        self.resample_labels()
        return self._subsampled_labeled_indices

    @subsampled_labeled_indices.setter
    def subsampled_labeled_indices(self, value: List[int]):
        """Sets the subsampled labeled indices."""
        self._subsampled_labeled_indices = value

    def resample_labels(self):
        """Resamples the labels."""
        subsampled = list(self.labeled_indices.values())
        if self.n_samples_per_label is not None:
            _subsampled = []
            for indices in subsampled:
                replace = len(indices) < self.n_samples_per_label
                subset = self.generator.choice(
                    indices, size=self.n_samples_per_label, replace=replace
                )
                _subsampled.append(subset)
            subsampled = _subsampled

        self.subsampled_labeled_indices = np.concatenate(subsampled)
        # triggers reassignment of dataloaders
        self.indices_list = [self.indices, self.subsampled_labeled_indices]
