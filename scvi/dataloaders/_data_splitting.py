import logging
from math import ceil, floor
from typing import Dict, List, Optional, Union

import lightning.pytorch as pl
import numpy as np
import torch
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
from scvi.utils._exceptions import InvalidParameterError

logger = logging.getLogger(__name__)


def validate_data_split(
    n_obs: int,
    train_size: Optional[float],
    validation_size: Optional[float],
    train_indices: Optional[List[int]],
    validation_indices: Optional[List[int]],
):
    """Validate data splitting parameters.

    Parameters
    ----------
    n_obs
        Number of observations in the dataset.
    train_size
        Fraction of the dataset for the training set.
    validation_size
        Fraction of the dataset for the validation set. If `None`, the validation set
        will be of size `1 - train_size`.
    train_indices
        Indices of the training set. Cannot be set if `train_size` is not `None`.
    validation_indices
        Indices of the validation set. Cannot be set if `validation_size` is not `None`.
        Set to the remaining indices if `train_indices` is not `None`. If the union of
        the two sets is not the full set of indices, the remaining indices are used for
        the test set.

    Returns
    -------
    The size of the training, validation, and test sets as a tuple if `train_size` is
    not `None`, otherwise the indices of the training, validation, and test sets as a
    tuple if `train_indices` is not `None`.
    """
    if train_size is None and train_indices is None:
        raise ValueError("Either `train_size` or `train_indices` must be specified.")
    if train_size is not None and train_indices is not None:
        raise ValueError("`train_size` and `train_indices` cannot both be specified.")

    if train_size is None and validation_size is not None:
        raise ValueError(
            "`train_size` must be specified if `validation_size` is specified."
        )
    if train_indices is None and validation_indices is not None:
        raise ValueError(
            "`train_indices` must be specified if `validation_indices` is specified."
        )

    if train_size is not None:
        if train_size > 1.0 or train_size <= 0.0:
            raise InvalidParameterError(
                "train_size",
                train_size,
                additional_message="`train_size` must be between 0 and 1.",
            )

        n_train = ceil(train_size * n_obs)

        if validation_size is None:
            n_val = n_obs - n_train
        elif validation_size >= 1.0 or validation_size < 0.0:
            raise InvalidParameterError(
                "validation_size",
                validation_size,
                additional_message="`validation_size` must be between 0 and 1.",
            )
        elif (train_size + validation_size) > 1:
            raise InvalidParameterError(
                "train_size + validation_size",
                train_size + validation_size,
                additional_message="`train_size + validation_size` must be between 0 and 1.",
            )
        else:
            n_val = floor(n_obs * validation_size)

        n_test = n_obs - n_train - n_val

        logging.info(
            f"Using {n_train} observations for training, {n_val} for validation "
            f"and {n_test} for testing."
        )
        return n_train, n_val, n_test

    if train_indices is not None:
        train_indices = np.array(train_indices)
        if validation_indices is not None:
            validation_indices = np.array(validation_indices)

        if np.amax(train_indices) >= n_obs or np.amin(train_indices) < 0:
            raise InvalidParameterError(
                "train_indices",
                train_indices,
                additional_message="`train_indices` contains invalid indices.",
            )

        if validation_indices is None:
            validation_indices = np.setdiff1d(np.arange(n_obs), train_indices)
        elif np.amax(validation_indices) >= n_obs or np.amin(validation_indices) < 0:
            raise InvalidParameterError(
                "validation_indices",
                validation_indices,
                additional_message="`validation_indices` contains invalid indices.",
            )

        union_indices = np.union1d(train_indices, validation_indices)
        test_indices = np.setdiff1d(np.arange(n_obs), union_indices)

        logging.info(
            f"Using {len(train_indices)} observations for training, "
            f"{len(validation_indices)} for validation and {len(test_indices)} for "
            "testing."
        )

        return train_indices, validation_indices, test_indices


class DataSplitter(pl.LightningDataModule):
    """Creates :class:`~scvi.data.AnnDataLoader` objects for train/validation/test sets.

    The test split is only created if `train_size + validation_size  < 1` or if the
    union of `train_indices` and `validation_indices` is not the full set of indices.

    Parameters
    ----------
    adata_manager
        :class:`~scvi.data.AnnDataManager` object that has been created via ``setup_anndata``.
    train_size
        Fraction of the dataset for the training set.
    validation_size
        Fraction of the dataset for the validation set. If `None`, the validation set
        will be of size `1 - train_size`.
    train_indices
        Indices of the training set. Ignored if `train_size` is not `None`.
    validation_indices
        Indices of the validation set. Ignored if `validation_size` is not `None`. Set
        to the remaining indices if `train_indices` is not `None`. If the union of the
        two sets is not the full set of indices, the remaining indices are used for the
        test set.
    shuffle
        Whether or not to shuffle the data before splitting. Ignored if `train_indices`
        is not `None`.
    pin_memory
        Whether to copy tensors into device-pinned memory before returning them. Passed
        into :class:`~scvi.data.AnnDataLoader`.
    **kwargs
        Keyword args for data loader. If adata has labeled data, data loader
        class is :class:`~scvi.dataloaders.SemiSupervisedDataLoader`,
        else data loader class is :class:`~scvi.dataloaders.AnnDataLoader`.

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
        if self.train_indices is not None:
            return

        all_indices = np.arange(self.adata_manager.adata.n_obs)
        if self.shuffle:
            random_state = np.random.default_rng(seed=settings.seed)
            all_indices = random_state.permutation(all_indices)

        self._train_dataloader = AnnDataLoader(
            self.adata_manager,
            indices=all_indices[: self.n_train],
            shuffle=True,
            drop_last=False,
            pin_memory=self.pin_memory,
            **self.data_loader_kwargs,
        )

        n_val_train = self.n_train + self.n_validation
        self._validation_dataloader = None
        self._test_dataloader = None

        if self.n_validation > 0:
            self._validation_dataloader = AnnDataLoader(
                self.adata_manager,
                indices=all_indices[self.n_train : n_val_train],
                shuffle=False,
                drop_last=False,
                pin_memory=self.pin_memory,
                **self.data_loader_kwargs,
            )
        if self.n_test > 0:
            self._test_dataloader = AnnDataLoader(
                self.adata_manager,
                indices=all_indices[n_val_train:],
                shuffle=False,
                drop_last=False,
                pin_memory=self.pin_memory,
                **self.data_loader_kwargs,
            )

    @property
    def train_dataloader(self):
        """Returns the train split data loader."""
        return self._train_dataloader

    @property
    def validation_dataloader(self):
        """Create validation split data loader."""
        return self._validation_dataloader

    @property
    def test_dataloader(self):
        """Create test split data loader."""
        return self._test_dataloader


class SemiSupervisedDataSplitter(pl.LightningDataModule):
    """Creates data loaders ``train_set``, ``validation_set``, ``test_set``.

    If ``train_size + validation_set < 1`` then ``test_set`` is non-empty.
    The ratio between labeled and unlabeled data in adata will be preserved
    in the train/test/val sets.

    Parameters
    ----------
    adata_manager
        :class:`~scvi.data.AnnDataManager` object that has been created via ``setup_anndata``.
    train_size
        float, or None (default is 0.9)
    validation_size
        float, or None (default is None)
    n_samples_per_label
        Number of subsamples for each label class to sample per epoch
    pin_memory
        Whether to copy tensors into device-pinned memory before returning them. Passed
        into :class:`~scvi.data.AnnDataLoader`.
    **kwargs
        Keyword args for data loader. If adata has labeled data, data loader
        class is :class:`~scvi.dataloaders.SemiSupervisedDataLoader`,
        else data loader class is :class:`~scvi.dataloaders.AnnDataLoader`.

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
        train_size: float = 0.9,
        validation_size: Optional[float] = None,
        n_samples_per_label: Optional[int] = None,
        pin_memory: bool = False,
        **kwargs,
    ):
        super().__init__()
        self.adata_manager = adata_manager
        self.train_size = float(train_size)
        self.validation_size = validation_size
        self.data_loader_kwargs = kwargs
        self.n_samples_per_label = n_samples_per_label

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

        self.data_loader_kwargs = kwargs
        self.pin_memory = pin_memory or settings.dl_pin_memory_gpu_training

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
