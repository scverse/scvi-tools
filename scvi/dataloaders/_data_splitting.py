from math import ceil, floor
from typing import Optional

import numpy as np
from anndata import AnnData

from scvi import _CONSTANTS, settings
from scvi.dataloaders._ann_dataloader import AnnDataLoader
from scvi.dataloaders._semi_dataloader import SemiSupervisedDataLoader
from scvi.model._utils import parse_use_gpu_arg


def validate_data_split(
    n_samples: int, train_size: float, validation_size: Optional[float] = None
):
    """
    Check data splitting parameters and return n_train and n_val.

    Parameters
    ----------
    n_samples
        Number of samples to split
    train_size
        Size of train set. Need to be: 0 < train_size <= 1.
    validation_size
        Size of validation set. Need to be 0 <= validation_size < 1
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
            "With n_samples={}, train_size={} and validation_size={}, the "
            "resulting train set will be empty. Adjust any of the "
            "aforementioned parameters.".format(n_samples, train_size, validation_size)
        )

    return n_train, n_val


class DataSplitter:
    """
    Creates data loaders ``train_set``, ``validation_set``, ``test_set``.

    If ``train_size + validation_set < 1`` then ``test_set`` is non-empty.

    Parameters
    ----------
    adata
        AnnData to split into train/test/val sets
    train_size
        float, or None (default is 0.9)
    validation_size
        float, or None (default is None)
    use_gpu
        Which GPU to use
    **kwargs
        Keyword args for data loader. If adata has labeled data, data loader
        class is :class:`~scvi.dataloaders.SemiSupervisedDataLoader`,
        else data loader class is :class:`~scvi.dataloaders.AnnDataLoader`.

    Examples
    --------
    >>> adata = scvi.data.synthetic_iid()
    >>> splitter = DataSplitter(adata)
    >>> train_dl, val_dl, test_dl = splitter()
    """

    def __init__(
        self,
        adata: AnnData,
        train_size: float = 0.9,
        validation_size: Optional[float] = None,
        use_gpu: bool = False,
        **kwargs,
    ):
        self.adata = adata
        self.train_size = float(train_size)
        self.validation_size = validation_size
        self.data_loader_kwargs = kwargs
        self.use_gpu = use_gpu

        self.train_idx, self.test_idx, self.val_idx = self.make_splits()

    def make_splits(self):
        """Split indices in train/test/val sets."""
        n = self.adata.n_obs
        n_train, n_val = validate_data_split(n, self.train_size, self.validation_size)
        random_state = np.random.RandomState(seed=settings.seed)
        permutation = random_state.permutation(n)
        val_idx = permutation[:n_val]
        train_idx = permutation[n_val : (n_val + n_train)]
        test_idx = permutation[(n_val + n_train) :]
        return train_idx, test_idx, val_idx

    def __call__(self, remake_splits=False):
        if remake_splits:
            self.train_idx, self.test_idx, self.val_idx = self.make_splits()

        gpus = parse_use_gpu_arg(self.use_gpu, return_device=False)
        pin_memory = (
            True if (settings.dl_pin_memory_gpu_training and gpus != 0) else False
        )

        # do not remove drop_last=3, skips over small minibatches
        return (
            AnnDataLoader(
                self.adata,
                indices=self.train_idx,
                shuffle=True,
                drop_last=3,
                pin_memory=pin_memory,
                **self.data_loader_kwargs,
            ),
            AnnDataLoader(
                self.adata,
                indices=self.val_idx,
                shuffle=True,
                drop_last=3,
                pin_memory=pin_memory,
                **self.data_loader_kwargs,
            ),
            AnnDataLoader(
                self.adata,
                indices=self.test_idx,
                shuffle=True,
                drop_last=3,
                pin_memory=pin_memory,
                **self.data_loader_kwargs,
            ),
        )


class SemiSupervisedDataSplitter:
    """
    Creates data loaders ``train_set``, ``validation_set``, ``test_set``.

    If ``train_size + validation_set < 1`` then ``test_set`` is non-empty.
    The ratio between labeled and unlabeled data in adata will be preserved
    in the train/test/val sets.

    Parameters
    ----------
    adata
        AnnData to split into train/test/val sets
    unlabeled_category
        Category to treat as unlabeled
    train_size
        float, or None (default is 0.9)
    validation_size
        float, or None (default is None)
    n_samples_per_label
        Number of subsamples for each label class to sample per epoch
    **kwargs
        Keyword args for data loader. If adata has labeled data, data loader
        class is :class:`~scvi.dataloaders.SemiSupervisedDataLoader`,
        else data loader class is :class:`~scvi.dataloaders.AnnDataLoader`.

    Examples
    --------
    >>> adata = scvi.data.synthetic_iid()
    >>> unknown_label = 'label_0'
    >>> splitter = SemiSupervisedDataSplitter(adata, unknown_label)
    >>> train_dl, val_dl, test_dl = splitter()
    """

    def __init__(
        self,
        adata: AnnData,
        unlabeled_category,
        train_size: float = 0.9,
        validation_size: Optional[float] = None,
        n_samples_per_label: Optional[int] = None,
        use_gpu: bool = False,
        **kwargs,
    ):
        self.adata = adata
        self.unlabeled_category = unlabeled_category
        self.train_size = float(train_size)
        self.validation_size = validation_size
        self.data_loader_kwargs = kwargs
        self.n_samples_per_label = n_samples_per_label

        setup_dict = adata.uns["_scvi"]
        key = setup_dict["data_registry"][_CONSTANTS.LABELS_KEY]["attr_key"]
        original_key = setup_dict["categorical_mappings"][key]["original_key"]
        labels = np.asarray(adata.obs[original_key]).ravel()
        self._unlabeled_indices = np.argwhere(labels == unlabeled_category).ravel()
        self._labeled_indices = np.argwhere(labels != unlabeled_category).ravel()

        self.data_loader_kwargs = kwargs
        self.use_gpu = use_gpu
        self.train_idx, self.test_idx, self.val_idx = self.make_splits()

    def make_splits(self):
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

        indices_train = indices_train.astype(int)
        indices_val = indices_val.astype(int)
        indices_test = indices_test.astype(int)
        return indices_train, indices_test, indices_val

    def __call__(self, remake_splits=False):
        if remake_splits:
            self.train_idx, self.test_idx, self.val_idx = self.make_splits()

        gpus = parse_use_gpu_arg(self.use_gpu, return_device=False)
        pin_memory = (
            True if (settings.dl_pin_memory_gpu_training and gpus != 0) else False
        )

        if len(self._labeled_indices) != 0:
            data_loader_class = SemiSupervisedDataLoader
            dl_kwargs = {
                "unlabeled_category": self.unlabeled_category,
                "n_samples_per_label": self.n_samples_per_label,
            }
        else:
            data_loader_class = AnnDataLoader
            dl_kwargs = {}

        dl_kwargs.update(self.data_loader_kwargs)

        scanvi_train_dl = data_loader_class(
            self.adata,
            indices=self.train_idx,
            shuffle=True,
            drop_last=3,
            pin_memory=pin_memory,
            **dl_kwargs,
        )
        scanvi_val_dl = data_loader_class(
            self.adata,
            indices=self.val_idx,
            shuffle=True,
            drop_last=3,
            pin_memory=pin_memory,
            **dl_kwargs,
        )
        scanvi_test_dl = data_loader_class(
            self.adata,
            indices=self.test_idx,
            shuffle=True,
            drop_last=3,
            pin_memory=pin_memory,
            **dl_kwargs,
        )

        return scanvi_train_dl, scanvi_val_dl, scanvi_test_dl
