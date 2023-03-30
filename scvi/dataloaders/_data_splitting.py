import logging
from typing import Dict, List, Optional, Union

import lightning.pytorch as pl
import numpy as np
import torch
from _docstrings import data_splitting_dsp
from _utils import validate_data_split
from torch.utils.data import (
    BatchSampler,
    DataLoader,
    Dataset,
    RandomSampler,
    SequentialSampler,
)

from scvi import REGISTRY_KEYS, settings
from scvi.data import AnnDataManager
from scvi.data._utils import get_anndata_attribute
from scvi.dataloaders._ann_dataloader import AnnDataLoader
from scvi.dataloaders._semi_dataloader import SemiSupervisedDataLoader
from scvi.model._utils import parse_device_args
from scvi.utils._docstrings import devices_dsp

logger = logging.getLogger(__name__)


@data_splitting_dsp.dedent
class DataSplitter(pl.LightningDataModule):
    """Creates :class:`~scvi.data.AnnDataLoader` objects for train/validation/test sets.

    Parameters
    ----------
    %(adata_manager)s
    %(train_size)s
    %(validation_size)s
    %(train_indices)s
    %(validation_indices)s
    %(shuffle)s
    %(pin_memory)s
    **kwargs
        Keyword arguments passed into :class:`~scvi.data.AnnDataLoader`.

    Examples
    --------
    >>> adata = scvi.data.synthetic_iid()
    >>> scvi.model.SCVI.setup_anndata(adata)
    >>> adata_manager = scvi.model.SCVI(adata).adata_manager
    >>> splitter = DataSplitter(adata)
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
        pin_memory: bool = False,
        **kwargs,
    ):
        super().__init__()
        self.adata_manager = adata_manager

        self.train_size = None
        self.validation_size = None
        self.test_size = None
        self.train_indices = None
        self.validation_indices = None
        self.test_indices = None
        self.shuffle = shuffle

        self.data_loader_kwargs = kwargs
        self.pin_memory = pin_memory or settings.dl_pin_memory_gpu_training

        splits = validate_data_split(
            adata_manager.adata.n_obs,
            train_size,
            validation_size,
            train_indices,
            validation_indices,
        )

        if train_size is not None:
            self.n_train, self.n_validation, self.n_test = splits
        else:
            self.train_indices, self.validation_indices, self.test_indices = splits
            self.n_train, self.n_validation, self.n_test = (
                len(self.train_indices),
                len(self.validation_indices),
                len(self.test_indices),
            )

    def setup(self, stage: Optional[str] = None):
        """Assign indices to train/validation/test splits if necessary."""
        if self.train_indices is None:
            all_indices = np.arange(self.adata_manager.adata.n_obs)
            if self.shuffle:
                random_state = np.random.default_rng(seed=settings.seed)
                all_indices = random_state.permutation(all_indices)

            n_val_train = self.n_train + self.n_validation
            self.train_indices = all_indices[: self.n_train]
            self.validation_indices = all_indices[self.n_train : n_val_train]
            self.test_indices = all_indices[n_val_train:]

        self._train_dataloader = AnnDataLoader(
            self.adata_manager,
            indices=self.train_indices,
            shuffle=True,
            drop_last=False,
            pin_memory=self.pin_memory,
            **self.data_loader_kwargs,
        )

        self._validation_dataloader = None
        self._test_dataloader = None

        if self.n_validation > 0:
            self._validation_dataloader = AnnDataLoader(
                self.adata_manager,
                indices=self.validation_indices,
                shuffle=False,
                drop_last=False,
                pin_memory=self.pin_memory,
                **self.data_loader_kwargs,
            )
        if self.n_test > 0:
            self._test_dataloader = AnnDataLoader(
                self.adata_manager,
                indices=self.test_indices,
                shuffle=False,
                drop_last=False,
                pin_memory=self.pin_memory,
                **self.data_loader_kwargs,
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


class SemiSupervisedDataSplitter(pl.LightningDataModule):
    """Creates :class:`~scvi.data.AnnDataLoader` objects for train/validation/test sets.

    Preserves the ratio between labeled and unlabeled data between the splits.

    Parameters
    ----------
    %(adata_manager)s
    %(train_size)s
    %(validation_size)s
    %(train_indices)s
    %(validation_indices)s
    %(shuffle)s
    n_samples_per_label
        Number of subsamples for each label class to sample per epoch.
    %(pin_memory)s
    **kwargs
        Keyword arguments passed into :class:`~scvi.data.AnnDataLoader` if there is no
        labeled data or :class:`~scvi.data.SemiSupervisedDataLoader` if there is labeled
        data.

    Examples
    --------
    >>> adata = scvi.data.synthetic_iid()
    >>> scvi.model.SCVI.setup_anndata(adata, labels_key="labels")
    >>> adata_manager = scvi.model.SCVI(adata).adata_manager
    >>> unknown_label = 'label_0'
    >>> splitter = SemiSupervisedDataSplitter(adata, unknown_label)
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
        **kwargs,
    ):
        super().__init__()
        self.adata_manager = adata_manager

        self.train_size = None
        self.validation_size = None
        self.test_size = None
        self.train_indices = None
        self.validation_indices = None
        self.test_indices = None
        self.shuffle = shuffle
        self.n_samples_per_label = n_samples_per_label

        self.data_loader_kwargs = kwargs
        self.pin_memory = pin_memory or settings.dl_pin_memory_gpu_training

        labels_state_registry = adata_manager.get_state_registry(
            REGISTRY_KEYS.LABELS_KEY
        )
        labels = get_anndata_attribute(
            adata_manager.adata,
            adata_manager.data_registry.labels.attr_name,
            labels_state_registry.original_key,
        ).ravel()
        self.unlabeled_category = labels_state_registry.unlabeled_category
        self._unlabeled_indices = np.argwhere(labels == self.unlabeled_category).ravel()
        self._labeled_indices = np.argwhere(labels != self.unlabeled_category).ravel()

    def setup(self, stage: Optional[str] = None):
        """Split indices in train/test/val sets."""
        n_labeled_idx = len(self._labeled_indices)
        n_unlabeled_idx = len(self._unlabeled_indices)

        if n_labeled_idx != 0:
            n_labeled_train, n_labeled_val = validate_data_split(
                n_labeled_idx, self.train_size, self.validation_size
            )
            rs = np.random.RandomState(seed=settings.seed)
            labeled_permutation = rs.choice(
                self._labeled_indices, len(self._labeled_indices), replace=False
            )
            labeled_idx_val = labeled_permutation[:n_labeled_val]
            labeled_idx_train = labeled_permutation[
                n_labeled_val : (n_labeled_val + n_labeled_train)
            ]
            labeled_idx_test = labeled_permutation[(n_labeled_val + n_labeled_train) :]
        else:
            labeled_idx_test = []
            labeled_idx_train = []
            labeled_idx_val = []

        if n_unlabeled_idx != 0:
            n_unlabeled_train, n_unlabeled_val = validate_data_split(
                n_unlabeled_idx, self.train_size, self.validation_size
            )
            rs = np.random.RandomState(seed=settings.seed)
            unlabeled_permutation = rs.choice(
                self._unlabeled_indices, len(self._unlabeled_indices)
            )
            unlabeled_idx_val = unlabeled_permutation[:n_unlabeled_val]
            unlabeled_idx_train = unlabeled_permutation[
                n_unlabeled_val : (n_unlabeled_val + n_unlabeled_train)
            ]
            unlabeled_idx_test = unlabeled_permutation[
                (n_unlabeled_val + n_unlabeled_train) :
            ]
        else:
            unlabeled_idx_train = []
            unlabeled_idx_val = []
            unlabeled_idx_test = []

        indices_train = np.concatenate((labeled_idx_train, unlabeled_idx_train))
        indices_val = np.concatenate((labeled_idx_val, unlabeled_idx_val))
        indices_test = np.concatenate((labeled_idx_test, unlabeled_idx_test))

        self.train_idx = indices_train.astype(int)
        self.val_idx = indices_val.astype(int)
        self.test_idx = indices_test.astype(int)

        if len(self._labeled_indices) != 0:
            self.data_loader_class = SemiSupervisedDataLoader
            dl_kwargs = {
                "n_samples_per_label": self.n_samples_per_label,
            }
        else:
            self.data_loader_class = AnnDataLoader
            dl_kwargs = {}

        self.data_loader_kwargs.update(dl_kwargs)

    def train_dataloader(self):
        """Create the train data loader."""
        return self.data_loader_class(
            self.adata_manager,
            indices=self.train_idx,
            shuffle=True,
            drop_last=False,
            pin_memory=self.pin_memory,
            **self.data_loader_kwargs,
        )

    def val_dataloader(self):
        """Create the validation data loader."""
        if len(self.val_idx) > 0:
            return self.data_loader_class(
                self.adata_manager,
                indices=self.val_idx,
                shuffle=False,
                drop_last=False,
                pin_memory=self.pin_memory,
                **self.data_loader_kwargs,
            )
        else:
            pass

    def test_dataloader(self):
        """Create the test data loader."""
        if len(self.test_idx) > 0:
            return self.data_loader_class(
                self.adata_manager,
                indices=self.test_idx,
                shuffle=False,
                drop_last=False,
                pin_memory=self.pin_memory,
                **self.data_loader_kwargs,
            )
        else:
            pass


@devices_dsp.dedent
class DeviceBackedDataSplitter(DataSplitter):
    """Creates loaders for data that is already on device, e.g., GPU.

    If ``train_size + validation_set < 1`` then ``test_set`` is non-empty.

    Parameters
    ----------
    adata_manager
        :class:`~scvi.data.AnnDataManager` object that has been created via ``setup_anndata``.
    train_size
        float, or None (default is 0.9)
    validation_size
        float, or None (default is None)
    %(param_accelerator)s
    %(param_device)s
    pin_memory
        Whether to copy tensors into device-pinned memory before returning them. Passed
        into :class:`~scvi.data.AnnDataLoader`.
    shuffle
        if ``True``, shuffles indices before sampling for training set
    shuffle_test_val
        Shuffle test and validation indices.
    batch_size
        batch size of each iteration. If `None`, do not minibatch

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
        train_size: float = 1.0,
        validation_size: Optional[float] = None,
        accelerator: str = "auto",
        device: Union[int, str] = "auto",
        pin_memory: bool = False,
        shuffle: bool = False,
        shuffle_test_val: bool = False,
        batch_size: Optional[int] = None,
        **kwargs,
    ):
        super().__init__(
            adata_manager=adata_manager,
            train_size=train_size,
            validation_size=validation_size,
            pin_memory=pin_memory,
            **kwargs,
        )
        self.batch_size = batch_size
        self.shuffle = shuffle
        self.shuffle_test_val = shuffle_test_val
        _, _, self.device = parse_device_args(
            accelerator=accelerator, devices=device, return_device="torch"
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
        dataset = _DeviceBackedDataset(tensor_dict)
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


class _DeviceBackedDataset(Dataset):
    def __init__(self, tensor_dict: Dict[str, torch.Tensor]):
        self.data = tensor_dict

    def __getitem__(self, idx: List[int]) -> Dict[str, torch.Tensor]:
        return_dict = {}
        for key, value in self.data.items():
            return_dict[key] = value[idx]

        return return_dict

    def __len__(self):
        for _, value in self.data.items():
            return len(value)
