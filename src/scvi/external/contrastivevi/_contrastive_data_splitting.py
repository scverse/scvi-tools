import numpy as np

from scvi import settings
from scvi.data import AnnDataManager
from scvi.dataloaders import DataSplitter
from scvi.dataloaders._data_splitting import (
    validate_data_split,
    validate_data_split_with_external_indexing,
)

from ._contrastive_dataloader import ContrastiveDataLoader


class ContrastiveDataSplitter(DataSplitter):
    """Creates ContrastiveDataLoader for training, validation, and test set.

    If ``train_size + validation_set < 1`` then ``test_set`` is non-empty.

    Parameters
    ----------
    adata_manager
        :class:`~scvi.data.AnnDataManager` object that has been created via
        ``setup_anndata``.
    background_indices: Indices for background samples in adata.
    target_indices: Indices for target samples in adata.
    train_size
        float, or None (default is 0.9)
    validation_size
        float, or None (default is None)
    shuffle_set_split
        Whether to shuffle indices before splitting. If `False`, the val, train, and
        test set are split in the sequential order of the data according to
        `validation_size` and `train_size` percentages.
    load_sparse_tensor
        If `True`, loads sparse CSR or CSC arrays in the input dataset as sparse
        :class:`~torch.Tensor` with the same layout. Can lead to significant
        speedups in transferring data to GPUs, depending on the sparsity of the data.
        Passed into :class:`~scvi.data.AnnDataLoader`.
    pin_memory
        Whether to copy tensors into device-pinned memory before returning them. Passed
        into :class:`~scvi.data.AnnDataLoader`.
    external_indexing
        A list of data split indices in the order of training, validation, and test sets.
        Validation and test set are not required and can be left empty.
        Note that per group (train,valid,test) it will cover both the background and target indices
    **kwargs
        Keyword args for data loader. Data loader class is
        :class:`~scvi.dataloaders.AnnDataLoader`.
    """

    def __init__(
        self,
        adata_manager: AnnDataManager,
        background_indices: list[int],
        target_indices: list[int],
        train_size: float | None = None,
        validation_size: float | None = None,
        shuffle_set_split: bool = True,
        load_sparse_tensor: bool = False,
        pin_memory: bool = False,
        external_indexing: list[np.array, np.array, np.array] | None = None,
        **kwargs,
    ) -> None:
        super().__init__(
            adata_manager=adata_manager,
            train_size=train_size,
            validation_size=validation_size,
            shuffle_set_split=shuffle_set_split,
            load_sparse_tensor=load_sparse_tensor,
            pin_memory=pin_memory,
            **kwargs,
        )
        self.background_indices = background_indices
        self.target_indices = target_indices
        self.external_indexing = external_indexing

        self.n_background = len(background_indices)
        self.n_target = len(target_indices)
        if external_indexing is None:
            self.n_background_train, self.n_background_val = validate_data_split(
                self.n_background,
                self.train_size,
                self.validation_size,
                self.data_loader_kwargs.pop("batch_size", settings.batch_size),
                self.drop_last,
                self.train_size_is_none,
            )
            self.n_target_train, self.n_target_val = validate_data_split(
                self.n_target,
                self.train_size,
                self.validation_size,
                self.data_loader_kwargs.pop("batch_size", settings.batch_size),
                self.drop_last,
                self.train_size_is_none,
            )
        else:
            # we need to intersect the external indexing given with the bg/target indices
            self.background_train_idx, self.background_val_idx, self.background_test_idx = (
                np.intersect1d(self.external_indexing[n], self.background_indices)
                for n in range(3)
            )
            self.n_background_train, self.n_background_val = (
                validate_data_split_with_external_indexing(
                    self.n_background,
                    [self.background_train_idx, self.background_val_idx, self.background_test_idx],
                    self.data_loader_kwargs.pop("batch_size", settings.batch_size),
                    self.drop_last,
                )
            )
            self.background_train_idx, self.background_val_idx, self.background_test_idx = (
                self.background_train_idx.tolist(),
                self.background_val_idx.tolist(),
                self.background_test_idx.tolist(),
            )

            self.target_train_idx, self.target_val_idx, self.target_test_idx = (
                np.intersect1d(self.external_indexing[n], self.target_indices) for n in range(3)
            )
            self.n_target_train, self.n_target_val = validate_data_split_with_external_indexing(
                self.n_target,
                [self.target_train_idx, self.target_val_idx, self.target_test_idx],
                self.data_loader_kwargs.pop("batch_size", settings.batch_size),
                self.drop_last,
            )
            self.target_train_idx, self.target_val_idx, self.target_test_idx = (
                self.target_train_idx.tolist(),
                self.target_val_idx.tolist(),
                self.target_test_idx.tolist(),
            )

        self.n_train = self.n_background_train + self.n_target_train
        self.n_val = self.n_background_val + self.n_target_val

    def setup(self, stage: str | None = None):
        """Split background and target indices into train/val/test sets."""
        background_indices = self.background_indices
        n_background_train = self.n_background_train
        n_background_val = self.n_background_val

        target_indices = self.target_indices
        n_target_train = self.n_target_train
        n_target_val = self.n_target_val

        # Need to separate to the external and non-external cases of the target/bg indices
        if self.external_indexing is None:
            if self.shuffle_set_split:
                random_state = np.random.RandomState(seed=settings.seed)
                background_indices = random_state.permutation(background_indices).tolist()
                target_indices = random_state.permutation(target_indices).tolist()

            self.background_val_idx = background_indices[:n_background_val]
            self.background_train_idx = background_indices[
                n_background_val : (n_background_val + n_background_train)
            ]
            self.background_test_idx = background_indices[
                (n_background_val + n_background_train) :
            ]

            self.target_val_idx = target_indices[:n_target_val]
            self.target_train_idx = target_indices[n_target_val : (n_target_val + n_target_train)]
            self.target_test_idx = target_indices[(n_target_val + n_target_train) :]

        self.val_idx = self.background_val_idx + self.target_val_idx
        self.train_idx = self.background_train_idx + self.target_train_idx
        self.test_idx = self.background_test_idx + self.target_test_idx

    def train_dataloader(self) -> ContrastiveDataLoader:
        return ContrastiveDataLoader(
            adata_manager=self.adata_manager,
            background_indices=self.background_train_idx,
            target_indices=self.target_train_idx,
            shuffle=True,
            drop_last=self.drop_last,
            load_sparse_tensor=self.load_sparse_tensor,
            pin_memory=self.pin_memory,
            **self.data_loader_kwargs,
        )

    def val_dataloader(self) -> ContrastiveDataLoader:
        if len(self.background_val_idx) > 0 and len(self.target_val_idx) > 0:
            return ContrastiveDataLoader(
                adata_manager=self.adata_manager,
                background_indices=self.background_val_idx,
                target_indices=self.target_val_idx,
                shuffle=False,
                drop_last=False,
                load_sparse_tensor=self.load_sparse_tensor,
                pin_memory=self.pin_memory,
                **self.data_loader_kwargs,
            )
        else:
            pass

    def test_dataloader(self) -> ContrastiveDataLoader:
        if len(self.background_test_idx) > 0 and len(self.target_test_idx) > 0:
            return ContrastiveDataLoader(
                adata_manager=self.adata_manager,
                background_indices=self.background_test_idx,
                target_indices=self.target_test_idx,
                shuffle=False,
                drop_last=False,
                load_sparse_tensor=self.load_sparse_tensor,
                pin_memory=self.pin_memory,
                **self.data_loader_kwargs,
            )
        else:
            pass
