import warnings
from math import ceil, floor

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


def validate_data_split(
    n_samples: int,
    train_size: float,
    validation_size: float | None = None,
    batch_size: int | None = None,
    drop_last: bool | int = False,
    train_size_is_none: bool | int = True,
):
    """Check data splitting parameters and return n_train and n_val.

    Parameters
    ----------
    n_samples
        Number of samples to split
    train_size
        Size of train set. Need to be: 0 < train_size <= 1.
    validation_size
        Size of validation set. Need to be 0 <= validation_size < 1
    batch_size
        batch size of each iteration. If `None`, do not minibatch
    drop_last
        drops last non-full batch
    train_size_is_none
        Whether the user did not explicitly input train_size
    """
    if train_size > 1.0 or train_size <= 0.0:
        raise ValueError("Invalid train_size. Must be: 0 < train_size <= 1")

    n_train = ceil(train_size * n_samples)

    if validation_size is None:
        n_val = n_samples - n_train
    elif validation_size >= 1.0 or validation_size < 0.0:
        raise ValueError("Invalid validation_size. Must be 0 <= validation_size < 1")
    elif (train_size + validation_size) > 1:
        raise ValueError("train_size + validation_size must be between 0 and 1")
    else:
        n_val = floor(n_samples * validation_size)

    if n_train == 0:
        raise ValueError(
            f"With n_samples={n_samples}, train_size={train_size} and "
            f"validation_size={validation_size}, the resulting train set will be empty. Adjust "
            "any of the aforementioned parameters."
        )

    if batch_size is not None:
        num_of_cells = n_train % batch_size
        if (num_of_cells < 3 and num_of_cells > 0) and drop_last is False:
            if not train_size_is_none:
                warnings.warn(
                    f"Last batch will have a small size of {num_of_cells} "
                    f"samples. Consider changing settings.batch_size or batch_size in model.train "
                    f"from currently {batch_size} to avoid errors during model training.",
                    UserWarning,
                    stacklevel=settings.warnings_stacklevel,
                )
            else:
                n_train -= num_of_cells
                if n_val > 0:
                    n_val += num_of_cells
                    warnings.warn(
                        f"{num_of_cells} cells moved from training set to validation set."
                        f" if you want to avoid it please use train_size parameter during train.",
                        UserWarning,
                        stacklevel=settings.warnings_stacklevel,
                    )

    return n_train, n_val


def validate_data_split_with_external_indexing(
    n_samples: int,
    external_indexing: list[np.array, np.array, np.array] | None = None,
    batch_size: int | None = None,
    drop_last: bool | int = False,
):
    """Check data splitting parameters and return n_train and n_val.

    Parameters
    ----------
    n_samples
        Number of samples to split
    external_indexing
        A list of data split indices in the order of training, validation, and test sets.
        Validation and test set are not required and can be left empty.
    batch_size
        batch size of each iteration. If `None`, do not minibatch
    drop_last
        drops last non-full batch
    """
    if not isinstance(external_indexing, list):
        raise ValueError("External indexing is not of list type")

    # validate the structure of it
    # make sure 3 elements exists and impute with None if not
    if len(external_indexing) == 0:
        external_indexing = [None, None, None]
    if len(external_indexing) == 1:
        external_indexing.append(None)
        external_indexing.append(None)
    if len(external_indexing) == 2:
        external_indexing.append(None)
    # (we can assume not all lists are given by user and impute the rest with empty arrays)
    external_indexing[0], external_indexing[1], external_indexing[2] = (
        np.array([]) if external_indexing[n] is None else external_indexing[n] for n in range(3)
    )
    if not all(isinstance(elem, np.ndarray) for elem in external_indexing):
        raise ValueError("One of the given external indexing arrays is not a np.array")

    # From this point on we will use the unique elements only
    external_indexing_unique = [
        set(external_indexing[0]),
        set(external_indexing[1]),
        set(external_indexing[2]),
    ]

    # check for duplications per subset
    if len(external_indexing_unique[0]) < len(external_indexing[0]):
        raise Warning("There are duplicate indexing in train set")
    if len(external_indexing_unique[1]) < len(external_indexing[1]):
        raise Warning("There are duplicate indexing in valid set")
    if len(external_indexing_unique[2]) < len(external_indexing[2]):
        raise Warning("There are duplicate indexing in test set")

    # check for total number of indexes (overlapping or missing)
    if (
        len(external_indexing_unique[0])
        + len(external_indexing_unique[1])
        + len(external_indexing_unique[2])
    ) < n_samples:
        raise Warning("There are missing indices please fix or remove those lines")

    if len(external_indexing_unique[0].intersection(external_indexing_unique[1])) != 0:
        raise ValueError("There are overlapping indexing between train and valid sets")
    if len(external_indexing_unique[0].intersection(external_indexing_unique[2])) != 0:
        raise ValueError("There are overlapping indexing between train and test sets")
    if len(external_indexing_unique[2].intersection(external_indexing_unique[1])) != 0:
        raise ValueError("There are overlapping indexing between test and valid sets")

    n_train = len(external_indexing[0])
    n_val = len(external_indexing[1])

    if batch_size is not None:
        num_of_cells = n_train % batch_size
        if (num_of_cells < 3 and num_of_cells > 0) and drop_last is False:
            warnings.warn(
                f"Last batch will have a small size of {num_of_cells} "
                f"samples. Consider changing settings.batch_size or batch_size in model.train "
                f"from currently {batch_size} to avoid errors during model training "
                f"or change the given external indices accordingly.",
                UserWarning,
                stacklevel=settings.warnings_stacklevel,
            )

    return n_train, n_val


class DataSplitter(pl.LightningDataModule):
    """Creates data loaders ``train_set``, ``validation_set``, ``test_set``.

    If ``train_size + validation_set < 1`` then ``test_set`` is non-empty.

    Parameters
    ----------
    adata_manager
        :class:`~scvi.data.AnnDataManager` object that has been created via ``setup_anndata``.
    train_size
        float, or None (default is None, which is practicaly 0.9 and potentially adding small last
        batch to validation cells)
    validation_size
        float, or None (default is None)
    shuffle_set_split
        Whether to shuffle indices before splitting. If `False`, the val, train, and test set are
        split in the sequential order of the data according to `validation_size` and `train_size`
        percentages.
    load_sparse_tensor
        ``EXPERIMENTAL`` If `True`, loads sparse CSR or CSC arrays in the input dataset as sparse
        :class:`~torch.Tensor` with the same layout. Can lead to significant speedups in
        transferring data to GPUs, depending on the sparsity of the data.
    pin_memory
        Whether to copy tensors into device-pinned memory before returning them. Passed
        into :class:`~scvi.data.AnnDataLoader`.
    external_indexing
        A list of data split indices in the order of training, validation, and test sets.
        Validation and test set are not required and can be left empty.
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

    data_loader_cls = AnnDataLoader

    def __init__(
        self,
        adata_manager: AnnDataManager,
        train_size: float | None = None,
        validation_size: float | None = None,
        shuffle_set_split: bool = True,
        load_sparse_tensor: bool = False,
        pin_memory: bool = False,
        external_indexing: list[np.array, np.array, np.array] | None = None,
        **kwargs,
    ):
        super().__init__()
        self.adata_manager = adata_manager
        self.train_size_is_none = not bool(train_size)
        self.train_size = 0.9 if self.train_size_is_none else float(train_size)
        self.validation_size = validation_size
        self.shuffle_set_split = shuffle_set_split
        self.load_sparse_tensor = load_sparse_tensor
        self.drop_last = kwargs.pop("drop_last", False)
        self.data_loader_kwargs = kwargs
        self.pin_memory = pin_memory
        self.external_indexing = external_indexing

        if self.external_indexing is not None:
            self.n_train, self.n_val = validate_data_split_with_external_indexing(
                self.adata_manager.adata.n_obs,
                self.external_indexing,
                self.data_loader_kwargs.get("batch_size", settings.batch_size),
                self.drop_last,
            )
        else:
            self.n_train, self.n_val = validate_data_split(
                self.adata_manager.adata.n_obs,
                self.train_size,
                self.validation_size,
                self.data_loader_kwargs.get("batch_size", settings.batch_size),
                self.drop_last,
                self.train_size_is_none,
            )

    def setup(self, stage: str | None = None):
        """Split indices in train/test/val sets."""
        if self.external_indexing is not None:
            # The structure and its order are guaranteed at this stage
            # (can include missing indexes for some group)
            self.train_idx = self.external_indexing[0]
            self.val_idx = self.external_indexing[1]
            self.test_idx = self.external_indexing[2]
        else:
            # just like it used to be w/o external indexing
            n_train = self.n_train
            n_val = self.n_val
            indices = np.arange(self.adata_manager.adata.n_obs)

            if self.shuffle_set_split:
                random_state = np.random.RandomState(seed=settings.seed)
                indices = random_state.permutation(indices)

            self.val_idx = indices[:n_val]
            self.train_idx = indices[n_val : (n_val + n_train)]
            self.test_idx = indices[(n_val + n_train) :]

    def train_dataloader(self):
        """Create train data loader."""
        return self.data_loader_cls(
            self.adata_manager,
            indices=self.train_idx,
            shuffle=True,
            drop_last=self.drop_last,
            load_sparse_tensor=self.load_sparse_tensor,
            pin_memory=self.pin_memory,
            **self.data_loader_kwargs,
        )

    def val_dataloader(self):
        """Create validation data loader."""
        if len(self.val_idx) > 0:
            return self.data_loader_cls(
                self.adata_manager,
                indices=self.val_idx,
                shuffle=False,
                drop_last=False,
                load_sparse_tensor=self.load_sparse_tensor,
                pin_memory=self.pin_memory,
                **self.data_loader_kwargs,
            )
        else:
            pass

    def test_dataloader(self):
        """Create test data loader."""
        if len(self.test_idx) > 0:
            return self.data_loader_cls(
                self.adata_manager,
                indices=self.test_idx,
                shuffle=False,
                drop_last=False,
                load_sparse_tensor=self.load_sparse_tensor,
                pin_memory=self.pin_memory,
                **self.data_loader_kwargs,
            )
        else:
            pass

    def on_after_batch_transfer(self, batch, dataloader_idx):
        """Converts sparse tensors to dense if necessary."""
        if self.load_sparse_tensor:
            for key, val in batch.items():
                layout = val.layout if isinstance(val, torch.Tensor) else None
                if layout is torch.sparse_csr or layout is torch.sparse_csc:
                    batch[key] = val.to_dense()

        return batch


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
        float, or None (default is None, which is practicaly 0.9 and potentially adding small last
        batch to validation cells)
    validation_size
        float, or None (default is None)
    shuffle_set_split
        Whether to shuffle indices before splitting. If `False`, the val, train, and test set
        are split in the sequential order of the data according to `validation_size` and
        `train_size` percentages.
    n_samples_per_label
        Number of subsamples for each label class to sample per epoch
    pin_memory
        Whether to copy tensors into device-pinned memory before returning them. Passed
        into :class:`~scvi.data.AnnDataLoader`.
    external_indexing
        A list of data split indices in the order of training, validation, and test sets.
        Validation and test set are not required and can be left empty.
        Note that per group (train,valid,test) it will cover both the labeled and unlebeled parts
    **kwargs
        Keyword args for data loader. If adata has labeled data, data loader
        class is :class:`~scvi.dataloaders.SemiSupervisedDataLoader`,
        else data loader class is :class:`~scvi.dataloaders.AnnDataLoader`.

    Examples
    --------
    >>> adata = scvi.data.synthetic_iid()
    >>> scvi.model.SCVI.setup_anndata(adata, labels_key="labels")
    >>> adata_manager = scvi.model.SCVI(adata).adata_manager
    >>> unknown_label = "label_0"
    >>> splitter = SemiSupervisedDataSplitter(adata, unknown_label)
    >>> splitter.setup()
    >>> train_dl = splitter.train_dataloader()
    """

    def __init__(
        self,
        adata_manager: AnnDataManager | None = None,
        datamodule: pl.LightningDataModule | None = None,
        train_size: float | None = None,
        validation_size: float | None = None,
        shuffle_set_split: bool = True,
        n_samples_per_label: int | None = None,
        pin_memory: bool = False,
        external_indexing: list[np.array, np.array, np.array] | None = None,
        **kwargs,
    ):
        super().__init__()
        self.adata_manager = adata_manager
        self.train_size_is_none = not bool(train_size)
        self.train_size = 0.9 if self.train_size_is_none else float(train_size)
        self.validation_size = validation_size
        self.shuffle_set_split = shuffle_set_split
        self.drop_last = kwargs.pop("drop_last", False)
        self.data_loader_kwargs = kwargs
        self.n_samples_per_label = n_samples_per_label

        labels_state_registry = adata_manager.get_state_registry(REGISTRY_KEYS.LABELS_KEY)
        labels = get_anndata_attribute(
            adata_manager.adata,
            adata_manager.data_registry.labels.attr_name,
            labels_state_registry.original_key,
            mod_key=getattr(self.adata_manager.data_registry.labels, "mod_key", None),
        ).ravel()
        self.unlabeled_category = labels_state_registry.unlabeled_category
        self._unlabeled_indices = np.argwhere(labels == self.unlabeled_category).ravel()
        self._labeled_indices = np.argwhere(labels != self.unlabeled_category).ravel()

        self.pin_memory = pin_memory
        self.external_indexing = external_indexing

    def setup(self, stage: str | None = None):
        """Split indices in train/test/val sets."""
        n_labeled_idx = len(self._labeled_indices)
        n_unlabeled_idx = len(self._unlabeled_indices)

        if n_labeled_idx != 0:
            # Need to separate to the external and non-external cases of the labeled indices
            if self.external_indexing is not None:
                # first we need to intersect the external indexing given with the labeled indices
                labeled_idx_train, labeled_idx_val, labeled_idx_test = (
                    np.intersect1d(self.external_indexing[n], self._labeled_indices)
                    for n in range(3)
                )
                n_labeled_train, n_labeled_val = validate_data_split_with_external_indexing(
                    n_labeled_idx,
                    [labeled_idx_train, labeled_idx_val, labeled_idx_test],
                    self.data_loader_kwargs.get("batch_size", settings.batch_size),
                    self.drop_last,
                )
            else:
                n_labeled_train, n_labeled_val = validate_data_split(
                    n_labeled_idx,
                    self.train_size,
                    self.validation_size,
                    self.data_loader_kwargs.get("batch_size", settings.batch_size),
                    self.drop_last,
                    self.train_size_is_none,
                )

                labeled_permutation = self._labeled_indices
                if self.shuffle_set_split:
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
            # Need to separate to the external and non-external cases of the unlabeled indices
            if self.external_indexing is not None:
                # we need to intersect the external indexing given with the labeled indices
                unlabeled_idx_train, unlabeled_idx_val, unlabeled_idx_test = (
                    np.intersect1d(self.external_indexing[n], self._unlabeled_indices)
                    for n in range(3)
                )
                n_unlabeled_train, n_unlabeled_val = validate_data_split_with_external_indexing(
                    n_unlabeled_idx,
                    [unlabeled_idx_train, unlabeled_idx_val, unlabeled_idx_test],
                    self.data_loader_kwargs.get("batch_size", settings.batch_size),
                    self.drop_last,
                )
            else:
                n_unlabeled_train, n_unlabeled_val = validate_data_split(
                    n_unlabeled_idx,
                    self.train_size,
                    self.validation_size,
                    self.data_loader_kwargs.get("batch_size", settings.batch_size),
                    self.drop_last,
                    self.train_size_is_none,
                )

                unlabeled_permutation = self._unlabeled_indices
                if self.shuffle_set_split:
                    rs = np.random.RandomState(seed=settings.seed)
                    unlabeled_permutation = rs.choice(
                        self._unlabeled_indices,
                        len(self._unlabeled_indices),
                        replace=False,
                    )

                unlabeled_idx_val = unlabeled_permutation[:n_unlabeled_val]
                unlabeled_idx_train = unlabeled_permutation[
                    n_unlabeled_val : (n_unlabeled_val + n_unlabeled_train)
                ]
                unlabeled_idx_test = unlabeled_permutation[(n_unlabeled_val + n_unlabeled_train) :]
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
            drop_last=self.drop_last,
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
        float, or None (default is None, which is practicaly 0.9 and potentially adding small last
        batch to validation cells)
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
        train_size: float | None = None,
        validation_size: float | None = None,
        accelerator: str = "auto",
        device: int | str = "auto",
        pin_memory: bool = False,
        shuffle: bool = False,
        shuffle_test_val: bool = False,
        batch_size: int | None = None,
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

    def setup(self, stage: str | None = None):
        """Create the train, validation, and test indices."""
        super().setup()

        if self.shuffle is False:
            self.train_idx = np.sort(self.train_idx)
            self.val_idx = np.sort(self.val_idx) if len(self.val_idx) > 0 else self.val_idx
            self.test_idx = np.sort(self.test_idx) if len(self.test_idx) > 0 else self.test_idx

        self.train_tensor_dict = self._get_tensor_dict(self.train_idx, device=self.device)
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

    def _make_dataloader(self, tensor_dict: dict[str, torch.Tensor], shuffle):
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
    def __init__(self, tensor_dict: dict[str, torch.Tensor]):
        self.data = tensor_dict

    def __getitem__(self, idx: list[int]) -> dict[str, torch.Tensor]:
        return_dict = {}
        for key, value in self.data.items():
            return_dict[key] = value[idx]

        return return_dict

    def __len__(self):
        for _, value in self.data.items():
            return len(value)
