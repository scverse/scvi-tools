import numpy as np

from scvi.data import AnnDataManager

from ._ann_dataloader import AnnDataLoader, labelled_indices_generator, subsample_labels
from ._concat_dataloader import ConcatDataLoader


class SemiSupervisedDataLoader(ConcatDataLoader):
    """DataLoader that supports semisupervised training.

    Parameters
    ----------
    adata_manager
        :class:`~scvi.data.AnnDataManager` object that has been created via ``setup_anndata``.
    n_samples_per_label
        Number of subsamples for each label class to sample per epoch. By default, there
        is no label subsampling.
    indices
        The indices of the observations in the adata to load
    shuffle
        Whether the data should be shuffled
    batch_size
        minibatch size to load each iteration
    data_and_attributes
        Dictionary with keys representing keys in data registry (`adata_manager.data_registry`)
        and value equal to desired numpy loading type (later made into torch tensor).
        If `None`, defaults to all registered data.
    data_loader_kwargs
        Keyword arguments for :class:`~torch.utils.data.DataLoader`
    """

    def __init__(
        self,
        adata_manager: AnnDataManager,
        n_samples_per_label: int | None = None,
        indices: list[int] | None = None,
        shuffle: bool = False,
        batch_size: int = 128,
        data_and_attributes: dict | None = None,
        drop_last: bool | int = False,
        **data_loader_kwargs,
    ):
        adata = adata_manager.adata
        if indices is None:
            indices = np.arange(adata.n_obs)

        self.indices = np.asarray(indices)

        if len(self.indices) == 0:
            return None

        self.n_samples_per_label = n_samples_per_label
        self.data_loader_kwargs = data_loader_kwargs

        self.labeled_locs, labelled_idx = labelled_indices_generator(
            adata_manager, indices, self.indices, self.n_samples_per_label
        )

        super().__init__(
            adata_manager=adata_manager,
            indices_list=[self.indices, labelled_idx],
            shuffle=shuffle,
            batch_size=batch_size,
            data_and_attributes=data_and_attributes,
            drop_last=drop_last,
            **self.data_loader_kwargs,
        )

    def resample_labels(self):
        """Resamples the labeled data."""
        labelled_idx = subsample_labels(self.labeled_locs, self.n_samples_per_label)
        # self.dataloaders[0] iterates over full_indices
        # self.dataloaders[1] iterates over the labelled_indices
        # change the indices of the labelled set
        self.dataloaders[1] = AnnDataLoader(
            self.adata_manager,
            indices=labelled_idx,
            shuffle=self._shuffle,
            batch_size=self._batch_size,
            data_and_attributes=self.data_and_attributes,
            drop_last=self._drop_last,
            **self.data_loader_kwargs,
        )
