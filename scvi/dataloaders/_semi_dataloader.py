from typing import List, Optional

import numpy as np
from anndata import AnnData

from scvi import _CONSTANTS

from ._concat_dataloader import ConcatDataLoader


class SemiSupervisedDataLoader(ConcatDataLoader):
    """
    DataLoader that supports semisupervised training.

    Parameters
    ----------
    adata
        AnnData object that have been registered via :func:`~scvi.data.setup_anndata`.
    unlabeled_category
        Category to treat as unlabeled
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
        Dictionary with keys representing keys in data registry (`adata.uns["_scvi"]`)
        and value equal to desired numpy loading type (later made into torch tensor).
        If `None`, defaults to all registered data.
    data_loader_kwargs
        Keyword arguments for :class:`~torch.utils.data.DataLoader`
    """

    def __init__(
        self,
        adata: AnnData,
        unlabeled_category: str,
        n_samples_per_label: Optional[int] = None,
        indices: Optional[List[int]] = None,
        shuffle: bool = False,
        batch_size: int = 128,
        data_and_attributes: Optional[dict] = None,
        **data_loader_kwargs,
    ):
        if indices is None:
            indices = np.arange(adata.n_obs)

        self.indices = indices

        if len(indices) == 0:
            return None

        self.n_samples_per_label = n_samples_per_label

        key = adata.uns["_scvi"]["data_registry"][_CONSTANTS.LABELS_KEY]["attr_key"]
        labels_obs_key = adata.uns["_scvi"]["categorical_mappings"][key]["original_key"]

        # save a nested list of the indices per labeled category
        self.labeled_locs = []
        labels = np.unique(adata.obs[labels_obs_key][indices])
        for label in labels:
            if label != unlabeled_category:
                label_loc_idx = np.where(adata.obs[labels_obs_key][indices] == label)[0]
                label_loc = indices[label_loc_idx]
                self.labeled_locs.append(label_loc)
        labelled_idx = self.subsample_labels()

        super().__init__(
            adata=adata,
            indices_list=[indices, labelled_idx],
            shuffle=shuffle,
            batch_size=batch_size,
            data_and_attributes=data_and_attributes,
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
