from typing import List, Optional, Union

import numpy as np

from scvi import REGISTRY_KEYS
from scvi.data import AnnDataManager
from scvi.data._utils import get_anndata_attribute

from ._concat_dataloader import ConcatDataLoader


class SemiSupervisedDataLoader(ConcatDataLoader):
    """
    DataLoader that supports semisupervised training.

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
        n_samples_per_label: Optional[int] = None,
        indices: Optional[List[int]] = None,
        shuffle: bool = False,
        batch_size: int = 128,
        data_and_attributes: Optional[dict] = None,
        drop_last: Union[bool, int] = False,
        **data_loader_kwargs,
    ):
        adata = adata_manager.adata
        if indices is None:
            indices = np.arange(adata.n_obs)

        self.indices = np.asarray(indices)

        if len(self.indices) == 0:
            return None

        self.n_samples_per_label = n_samples_per_label

        labels_state_registry = adata_manager.get_state_registry(
            REGISTRY_KEYS.LABELS_KEY
        )
        labels = get_anndata_attribute(
            adata_manager.adata,
            adata_manager.data_registry.labels.attr_name,
            labels_state_registry.original_key,
        ).ravel()

        # save a nested list of the indices per labeled category
        self.labeled_locs = []
        for label in np.unique(labels):
            if label != labels_state_registry.unlabeled_category:
                label_loc_idx = np.where(labels[indices] == label)[0]
                label_loc = self.indices[label_loc_idx]
                self.labeled_locs.append(label_loc)
        labelled_idx = self.subsample_labels()

        super().__init__(
            adata_manager=adata_manager,
            indices_list=[self.indices, labelled_idx],
            shuffle=shuffle,
            batch_size=batch_size,
            data_and_attributes=data_and_attributes,
            drop_last=drop_last,
            **data_loader_kwargs,
        )

    def resample_labels(self):
        """Resamples the labeled data."""
        labelled_idx = self.subsample_labels()
        # self.dataloaders[0] iterates over full_indices
        # self.dataloaders[1] iterates over the labelled_indices
        # change the indices of the labelled set
        self.dataloaders[1].indices = labelled_idx

    def subsample_labels(self):
        """Subsamples each label class by taking up to n_samples_per_label samples per class."""
        if self.n_samples_per_label is None:
            return np.concatenate(self.labeled_locs)

        sample_idx = []
        for loc in self.labeled_locs:
            if len(loc) < self.n_samples_per_label:
                sample_idx.append(loc)
            else:
                label_subset = np.random.choice(
                    loc, self.n_samples_per_label, replace=False
                )
                sample_idx.append(label_subset)
        sample_idx = np.concatenate(sample_idx)
        return sample_idx
