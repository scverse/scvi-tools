import logging
from typing import Dict, List, Optional, Union

import lightning.pytorch as pl
import numpy as np
import torch
from torch.utils.data import (
    BatchSampler,
    DataLoader,
    RandomSampler,
    SequentialSampler,
)

from scvi import REGISTRY_KEYS, settings
from scvi.data import AnnDataManager
from scvi.data._utils import get_anndata_attribute
from scvi.model._utils import parse_device_args

from ._dataloaders import AnnDataLoader, SemiSupervisedDataLoader
from ._datasets import DeviceBackedDataset
from ._docstrings import data_splitting_dsp
from ._utils import validate_data_split

logger = logging.getLogger(__name__)


@data_splitting_dsp.dedent
class DataSplitter(pl.LightningDataModule):
    """%(summary)s

    Parameters
    ----------
    %(param_adata_manager)s
    %(param_train_size)s
    %(param_validation_size)s
    %(param_train_indices)s
    %(param_validation_indices)s
    %(param_shuffle)s
    %(param_pin_memory)s
    %(param_train_dataloader_kwargs)s
    %(param_validation_dataloader_kwargs)s
    %(param_test_dataloader_kwargs)s

    Examples
    --------
    >>> adata = scvi.data.synthetic_iid()
    >>> scvi.model.SCVI.setup_anndata(adata)
    >>> adata_manager = scvi.model.SCVI(adata).adata_manager
    >>> splitter = DataSplitter(adata)
    >>> splitter.setup()
    >>> train_dl = splitter.train_dataloader()
    """

    _data_loader_cls = AnnDataLoader

    def __init__(
        self,
        adata_manager: AnnDataManager,
        train_size: Optional[float] = 0.9,
        validation_size: Optional[float] = None,
        train_indices: Optional[List[int]] = None,
        validation_indices: Optional[List[int]] = None,
        shuffle: bool = True,
        pin_memory: bool = False,
        train_dataloader_kwargs: Optional[Dict] = None,
        validation_dataloader_kwargs: Optional[Dict] = None,
        test_dataloader_kwargs: Optional[Dict] = None,
    ):
        super().__init__()
        self.adata_manager = adata_manager
        pin_memory = pin_memory or settings.dl_pin_memory_gpu_training

        self.train_dataloader_kwargs = {
            "shuffle": True,
            "drop_last": False,
            "pin_memory": pin_memory,
        }
        self.train_dataloader_kwargs.update(train_dataloader_kwargs or {})
        self.validation_dataloader_kwargs = {
            "shuffle": False,
            "drop_last": False,
            "pin_memory": pin_memory,
        }
        self.validation_dataloader_kwargs.update(validation_dataloader_kwargs or {})
        self.test_dataloader_kwargs = {
            "shuffle": False,
            "drop_last": False,
            "pin_memory": pin_memory,
        }
        self.test_dataloader_kwargs.update(test_dataloader_kwargs or {})

        self.indices = validate_data_split(
            n_obs=adata_manager.adata.n_obs,
            train_size=train_size,
            validation_size=validation_size,
            train_indices=train_indices,
            validation_indices=validation_indices,
            shuffle=shuffle,
        )

    def setup(self, stage: Optional[str] = None):
        """Assign indices to train/validation/test splits if necessary."""
        self._train_dataloader = self._data_loader_cls(
            self.adata_manager,
            indices=self.indices.train,
            **self.train_dataloader_kwargs,
        )
        self._validation_dataloader = self._data_loader_cls(
            self.adata_manager,
            indices=self.indices.validation,
            **self.validation_dataloader_kwargs,
        )
        self._test_dataloader = self._data_loader_cls(
            self.adata_manager,
            indices=self.indices.test,
            **self.test_dataloader_kwargs,
        )

    def train_dataloader(self):
        """Returns the train split data loader."""
        return self._train_dataloader

    def val_dataloader(self):
        """Returns the validation split data loader."""
        return self._validation_dataloader

    def test_dataloader(self):
        """Returns the test split data loader."""
        return self._test_dataloader


@data_splitting_dsp.dedent
class SemiSupervisedDataSplitter(DataSplitter):
    """%(summary)s

    Preserves the ratio between labeled and unlabeled data between the splits.

    Parameters
    ----------
    %(param_adata_manager)s
    %(param_train_size)s
    %(param_validation_size)s
    %(param_train_indices)s
    %(param_validation_indices)s
    %(param_shuffle)s
    %(n_samples_per_label)s
    %(pin_memory)s
    %(param_train_dataloader_kwargs)s
    %(param_validation_dataloader_kwargs)s
    %(param_test_dataloader_kwargs)s

    Examples
    --------
    >>> adata = scvi.data.synthetic_iid()
    >>> scvi.model.SCVI.setup_anndata(adata, labels_key="labels")
    >>> adata_manager = scvi.model.SCVI(adata).adata_manager
    >>> splitter = SemiSupervisedDataSplitter(adata_manager)
    >>> splitter.setup()
    >>> train_dl = splitter.train_dataloader()
    """

    def __init__(
        self,
        adata_manager: AnnDataManager,
        train_size: Optional[float] = 0.9,
        validation_size: Optional[float] = None,
        train_indices: Optional[List[int]] = None,
        validation_indices: Optional[List[int]] = None,
        shuffle: bool = True,
        n_samples_per_label: Optional[int] = None,
        pin_memory: bool = False,
        train_dataloader_kwargs: Optional[Dict] = None,
        validation_dataloader_kwargs: Optional[Dict] = None,
        test_dataloader_kwargs: Optional[Dict] = None,
    ):
        super().__init__(
            adata_manager,
            train_size=train_size,
            validation_size=validation_size,
            train_indices=train_indices,
            validation_indices=validation_indices,
            shuffle=shuffle,
            pin_memory=pin_memory,
            train_dataloader_kwargs=train_dataloader_kwargs,
            validation_dataloader_kwargs=validation_dataloader_kwargs,
            test_dataloader_kwargs=test_dataloader_kwargs,
        )

        labels_state_registry = self.adata_manager.get_state_registry(
            REGISTRY_KEYS.LABELS_KEY
        )
        labels = get_anndata_attribute(
            self.adata_manager.adata,
            self.adata_manager.data_registry.labels.attr_name,
            labels_state_registry.original_key,
        ).ravel()
        unlabeled_category = labels_state_registry.unlabeled_category
        unlabled_indices = np.argwhere(labels == unlabeled_category).ravel()
        labeled_indices = np.argwhere(labels != unlabeled_category).ravel()

        if len(labeled_indices) != 0:
            self.train_dataloader_kwargs["n_samples_per_label"] = n_samples_per_label
            self.validation_dataloader_kwargs[
                "n_samples_per_label"
            ] = n_samples_per_label
            self.test_dataloader_kwargs["n_samples_per_label"] = n_samples_per_label
            self._data_loader_cls = SemiSupervisedDataLoader

        labeled_indices = validate_data_split(
            all_indices=labeled_indices,
            train_size=train_size,
            validation_size=validation_size,
            train_indices=train_indices,
            validation_indices=validation_indices,
            shuffle=shuffle,
        )
        unlabeled_indices = validate_data_split(
            all_indices=unlabled_indices,
            train_size=train_size,
            validation_size=validation_size,
            train_indices=train_indices,
            validation_indices=validation_indices,
            shuffle=shuffle,
        )

        self.indices = labeled_indices + unlabeled_indices


@data_splitting_dsp
class DeviceBackedDataSplitter(DataSplitter):
    """%(summary)s

    Moves data to the device specified by `accelerator` and `device`.

    Parameters
    ----------
    %(param_adata_manager)s
    %(param_train_size)s
    %(param_validation_size)s
    %(param_train_indices)s
    %(param_validation_indices)s
    %(param_shuffle)s
    %(param_accelerator)s
    %(param_device)s
    %(param_pin_memory)s
    %(param_train_dataloader_kwargs)s
    %(param_validation_dataloader_kwargs)s
    %(param_test_dataloader_kwargs)s

    Examples
    --------
    >>> adata = scvi.data.synthetic_iid()
    >>> scvi.model.SCVI.setup_anndata(adata)
    >>> adata_manager = scvi.model.SCVI(adata).adata_manager
    >>> splitter = DeviceBackedDataSplitter(adata)
    >>> splitter.setup()
    >>> train_dl = splitter.train_dataloader()
    """

    def __init__(
        self,
        adata_manager: AnnDataManager,
        train_size: Optional[float] = 0.9,
        validation_size: Optional[float] = None,
        train_indices: Optional[List[int]] = None,
        validation_indices: Optional[List[int]] = None,
        shuffle: bool = True,
        accelerator: str = "auto",
        device: Union[int, str] = "auto",
        pin_memory: bool = False,
        train_dataloader_kwargs: Optional[Dict] = None,
        validation_dataloader_kwargs: Optional[Dict] = None,
        test_dataloader_kwargs: Optional[Dict] = None,
    ):
        super().__init__(
            adata_manager,
            train_size=train_size,
            validation_size=validation_size,
            train_indices=train_indices,
            validation_indices=validation_indices,
            shuffle=shuffle,
            pin_memory=pin_memory,
            train_dataloader_kwargs=train_dataloader_kwargs,
            validation_dataloader_kwargs=validation_dataloader_kwargs,
            test_dataloader_kwargs=test_dataloader_kwargs,
        )
        _, _, self.device = parse_device_args(
            accelerator=accelerator,
            devices=device,
            return_device="torch",
            validate_single_device=True,
        )

    def setup(self, stage: Optional[str] = None):
        """Create the train, validation, and test indices."""
        super().setup()

        if self.shuffle is False:
            self.train_idx = np.sort(self.train_idx)
            self.val_idx = (
                np.sort(self.val_idx) if len(self.val_idx) > 0 else self.val_idx
            )
            self.test_idx = (
                np.sort(self.test_idx) if len(self.test_idx) > 0 else self.test_idx
            )

        self.train_tensor_dict = self._get_tensor_dict(
            self.train_idx, device=self.device
        )
        self.test_tensor_dict = self._get_tensor_dict(self.test_idx, device=self.device)
        self.val_tensor_dict = self._get_tensor_dict(self.val_idx, device=self.device)

    def _get_tensor_dict(self, indices, device):
        """Get tensor dict for a given set of indices."""
        if len(indices) is not None and len(indices) > 0:
            dl = AnnDataLoader(
                self.adata_manager,
                indices=indices,
                batch_size=len(indices),
                shuffle=False,
                pin_memory=self.pin_memory,
                **self.data_loader_kwargs,
            )
            # will only have one minibatch
            for batch in dl:
                tensor_dict = batch

            for k, v in tensor_dict.items():
                tensor_dict[k] = v.to(device)

            return tensor_dict
        else:
            return None

    def _make_dataloader(self, tensor_dict: Dict[str, torch.Tensor], shuffle):
        """Create a dataloader from a tensor dict."""
        if tensor_dict is None:
            return None
        dataset = DeviceBackedDataset(tensor_dict)
        bs = self.batch_size if self.batch_size is not None else len(dataset)
        sampler_cls = SequentialSampler if not shuffle else RandomSampler
        sampler = BatchSampler(
            sampler=sampler_cls(dataset),
            batch_size=bs,
            drop_last=False,
        )
        return DataLoader(dataset, sampler=sampler, batch_size=None)

    def train_dataloader(self):
        """Create the train data loader."""
        return self._make_dataloader(self.train_tensor_dict, self.shuffle)

    def test_dataloader(self):
        """Create the test data loader."""
        return self._make_dataloader(self.test_tensor_dict, self.shuffle_test_val)

    def val_dataloader(self):
        """Create the validation data loader."""
        return self._make_dataloader(self.val_tensor_dict, self.shuffle_test_val)
