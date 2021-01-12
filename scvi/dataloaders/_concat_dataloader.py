from typing import Optional

import numpy as np
from torch.utils.data import DataLoader
from ._scvi_dataloader import ScviDataLoader
from itertools import cycle


class ConcatDataLoader(DataLoader):
    def __init__(
        self,
        adata,
        indices_list,
        shuffle=False,
        batch_size=128,
        data_and_attributes: Optional[dict] = None,
        **data_loader_kwargs,
    ):
        self.dataloaders = []
        for indices in indices_list:
            self.dataloaders.append(
                ScviDataLoader(
                    adata,
                    indices=indices,
                    shuffle=shuffle,
                    batch_size=batch_size,
                    data_and_attributes=data_and_attributes,
                    **data_loader_kwargs,
                )
            )
        lens = [len(dl) for dl in self.dataloaders]
        self.largest_dl = self.dataloaders[np.argmax(lens)]
        super().__init__(self.largest_dl, **data_loader_kwargs)

    def __len__(self):
        return len(self.largest_dl)

    def __iter__(self):
        iter_list = [
            cycle(dl) if dl != self.largest_dl else dl for dl in self.dataloaders
        ]
        return zip(*iter_list)


class SemiSupervisedDataLoader(ConcatDataLoader):
    def __init__(
        self,
        adata,
        labels_obs_key,
        unlabeled_category,
        n_samples_per_label,
        indices=None,
        shuffle=False,
        batch_size=128,
        data_and_attributes: Optional[dict] = None,
        **data_loader_kwargs,
    ):
        if indices is None:
            indices = np.arange(adata.n_obs)

        self.indices = indices

        if len(indices) == 0:
            return None

        self.n_samples_per_label = n_samples_per_label

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
        labelled_idx = self.subsample_labels()
        # self.dataloaders[0] iterates over full_indices
        # self.dataloaders[1] iterates over the labelled_indices
        # change the indicees of the labelled set
        self.dataloaders[1].indices = labelled_idx

    def subsample_labels(self):
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
