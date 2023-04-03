import logging
from typing import Dict, List, Optional, Union

import lightning.pytorch as pl
import numpy as np

from scvi import REGISTRY_KEYS, settings
from scvi.data import AnnDataManager
from scvi.data._utils import get_anndata_attribute
from scvi.model._utils import parse_device_args

from ._dataloaders import (
    AnnDataLoader,
    DeviceBackedDataLoader,
    SemiSupervisedDataLoader,
)
from ._docstrings import dataloaders_dsp
from ._utils import validate_data_split

logger = logging.getLogger(__name__)


@dataloaders_dsp.dedent
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


@dataloaders_dsp.dedent
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


@dataloaders_dsp.dedent
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

    _data_loader_cls = DeviceBackedDataLoader

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
        _, _, self.device = parse_device_args(
            accelerator=accelerator,
            devices=device,
            return_device="torch",
            validate_single_device=True,
        )
        super().__init__(
            adata_manager,
            train_size=train_size,
            validation_size=validation_size,
            train_indices=train_indices,
            validation_indices=validation_indices,
            shuffle=shuffle,
            pin_memory=pin_memory,
            pin_memory_device=self.device,
            train_dataloader_kwargs=train_dataloader_kwargs,
            validation_dataloader_kwargs=validation_dataloader_kwargs,
            test_dataloader_kwargs=test_dataloader_kwargs,
        )

    def setup(self, stage: Optional[str] = None):
        """Assign indices to train/validation/test splits if necessary."""
        self._train_dataloader = self._data_loader_cls(
            self.adata_manager,
            self.device,
            indices=self.indices.train,
            **self.train_dataloader_kwargs,
        )
        self._validation_dataloader = self._data_loader_cls(
            self.adata_manager,
            self.device,
            indices=self.indices.validation,
            **self.validation_dataloader_kwargs,
        )
        self._test_dataloader = self._data_loader_cls(
            self.adata_manager,
            self.device,
            indices=self.indices.test,
            **self.test_dataloader_kwargs,
        )
