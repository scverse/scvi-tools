from typing import Optional

import numpy as np
from anndata import AnnData
from sklearn.model_selection._split import _validate_shuffle_split

from scvi import _CONSTANTS, settings
from scvi.dataloaders._ann_dataloader import AnnDataLoader
from scvi.dataloaders._semi_dataloader import SemiSupervisedDataLoader


class DataSplitter:
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
        if self.train_size > 1.0 or self.train_size <= 0.0:
            raise ValueError(
                "train_size needs to be greater than 0 and less than or equal to 1"
            )
        try:
            n = self.adata.n_obs
            n_train, n_val = _validate_shuffle_split(
                n, self.validation_size, self.train_size
            )
        except ValueError:
            if self.train_size != 1.0:
                raise ValueError(
                    "Choice of train_size={} and validation_size={} not understood".format(
                        self.train_size, self.validation_size
                    )
                )
            n_train, n_val = n, 0
        random_state = np.random.RandomState(seed=settings.seed)
        permutation = random_state.permutation(n)
        val_idx = permutation[:n_val]
        train_idx = permutation[n_val : (n_val + n_train)]
        test_idx = permutation[(n_val + n_train) :]
        return train_idx, test_idx, val_idx

    def __call__(self, remake_splits=False):
        if remake_splits:
            self.train_idx, self.test_idx, self.val_idx = self.make_splits()

        pin_memory = (
            True if (settings.dl_pin_memory_gpu_training and self.use_gpu) else False
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
    """

    def __init__(
        self,
        adata: AnnData,
        unlabeled_category,
        train_size: float = 0.9,
        validation_size: Optional[float] = None,
        n_samples_per_label: Optional[int] = None,
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
        self.train_idx, self.test_idx, self.val_idx = self.make_splits()

    def make_splits(self):
        if self.train_size > 1.0 or self.train_size <= 0.0:
            raise ValueError(
                "train_size needs to be greater than 0 and less than or equal to 1"
            )

        n_labeled_idx = len(self._labeled_indices)
        n_unlabeled_idx = len(self._unlabeled_indices)

        def get_train_val_split(n_samples, test_size, train_size):
            try:
                n_train, n_val = _validate_shuffle_split(
                    n_samples, test_size, train_size
                )
            except ValueError:
                if train_size != 1.0 and n_samples != 1:
                    raise ValueError(
                        "Choice of train_size={} and validation_size={} not understood".format(
                            train_size, test_size
                        )
                    )
                n_train, n_val = n_samples, 0
            return n_train, n_val

        if n_labeled_idx != 0:
            n_labeled_train, n_labeled_val = get_train_val_split(
                n_labeled_idx, self.validation_size, self.train_size
            )
            labeled_permutation = np.random.choice(
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
            n_unlabeled_train, n_unlabeled_val = get_train_val_split(
                n_unlabeled_idx, self.validation_size, self.train_size
            )
            unlabeled_permutation = np.random.choice(
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
            **dl_kwargs,
        )
        scanvi_val_dl = data_loader_class(
            self.adata,
            indices=self.val_idx,
            shuffle=True,
            drop_last=3,
            **dl_kwargs,
        )
        scanvi_test_dl = data_loader_class(
            self.adata,
            indices=self.test_idx,
            shuffle=True,
            drop_last=3,
            **dl_kwargs,
        )

        return scanvi_train_dl, scanvi_val_dl, scanvi_test_dl
