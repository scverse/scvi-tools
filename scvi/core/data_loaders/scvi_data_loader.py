import copy
import logging
from typing import Optional
from scvi._compat import Literal
from itertools import cycle

import anndata
import numpy as np
import torch
from torch.utils.data import DataLoader

from scvi.data._scvidataset import ScviDataset

logger = logging.getLogger(__name__)


class BatchSampler(torch.utils.data.sampler.Sampler):
    """
    Custom torch Sampler that returns a list of indices of size batch_size.

    Parameters
    ----------
    indices
        list of indices to sample from
    batch_size
        batch size of each iteration
    shuffle
        if ``True``, shuffles indices before sampling

    """

    def __init__(
        self,
        indices: np.ndarray,
        batch_size: int,
        shuffle: bool,
        drop_last: bool = False,
    ):
        self.indices = indices
        self.n_obs = len(indices)
        self.batch_size = batch_size
        self.shuffle = shuffle

    def __iter__(self):
        if self.shuffle is True:
            idx = torch.randperm(len(self.indices)).tolist()
        else:
            idx = torch.arange(len(self.indices)).tolist()

        data_iter = iter(
            [
                self.indices[idx[i : i + self.batch_size]]
                for i in range(0, len(idx), self.batch_size)
            ]
        )
        return data_iter

    def __len__(self):
        from math import ceil

        length = ceil(self.n_obs / self.batch_size)

        return length


class ScviDataLoader(DataLoader):
    """
    Scvi Data Loader.

    A `ScviDataLoader` instance is instantiated with a model and a gene_dataset, and
    as well as additional arguments that for Pytorch's `DataLoader`. A subset of indices can be specified, for
    purposes such as splitting the data into train/test or labelled/unlabelled (for semi-supervised learning).
    Each trainer instance of the `Trainer` class can therefore have multiple `ScviDataLoader` instances to train a model.
    A `ScviDataLoader` instance also comes with methods to compute likelihood and other relevant training metrics.

    Parameters
    ----------
    adata
        An anndata instance
    shuffle
        Specifies if a `RandomSampler` or a `SequentialSampler` should be used
    indices
        Specifies how the data should be split with regards to train/test or labelled/unlabelled
    use_cuda
        Default: ``True``
    data_loader_kwargs
        Keyword arguments to passed into the `DataLoader`

    """

    def __init__(
        self,
        adata: anndata.AnnData,
        shuffle=False,
        indices=None,
        use_cuda=True,
        batch_size=128,
        data_and_attributes: Optional[dict] = None,
        **data_loader_kwargs,
    ):

        if "_scvi" not in adata.uns.keys():
            raise ValueError("Please run setup_anndata() on your anndata object first.")

        if data_and_attributes is not None:
            data_registry = adata.uns["_scvi"]["data_registry"]
            for key in data_and_attributes.keys():
                if key not in data_registry.keys():
                    raise ValueError(
                        "{} required for model but not included when setup_anndata was run".format(
                            key
                        )
                    )

        self.dataset = ScviDataset(adata, getitem_tensors=data_and_attributes)
        self.to_monitor = []
        self.use_cuda = use_cuda

        sampler_kwargs = {
            "batch_size": batch_size,
            "shuffle": shuffle,
        }

        if indices is None:
            indices = np.arange(len(self.dataset))
            sampler_kwargs["indices"] = indices
        else:
            if hasattr(indices, "dtype") and indices.dtype is np.dtype("bool"):
                indices = np.where(indices)[0].ravel()
            indices = np.asarray(indices)
            sampler_kwargs["indices"] = indices

        self.indices = indices
        self.sampler_kwargs = sampler_kwargs
        sampler = BatchSampler(**self.sampler_kwargs)
        self.data_loader_kwargs = copy.copy(data_loader_kwargs)
        # do not touch batch size here, sampler gives batched indices
        self.data_loader_kwargs.update({"sampler": sampler, "batch_size": None})

        super().__init__(self.dataset, **self.data_loader_kwargs)


class SemiSupervisedDataloader(ScviDataLoader):
    def __init__(
        self,
        adata,
        labeled_indices,
        unlabelled_indices,
        scheme: Literal["joint", "alternate", "both"] = "both",
        shuffle=False,
        use_cuda=True,
        batch_size=128,
        data_and_attributes: Optional[dict] = None,
        **data_loader_kwargs,
    ):
        self.full_dataset = ScviDataLoader(
            adata,
            shuffle=True,
            use_cuda=use_cuda,
            batch_size=batch_size,
            data_and_attributes=data_and_attributes,
            **data_loader_kwargs,
        )
        self.labelled_set = ScviDataLoader(
            adata,
            indices=labeled_indices,
            shuffle=shuffle,
            use_cuda=use_cuda,
            batch_size=batch_size,
            data_and_attributes=data_and_attributes,
            **data_loader_kwargs,
        )
        self.unlabelled_set = ScviDataLoader(
            adata,
            indices=labeled_indices,
            shuffle=shuffle,
            use_cuda=use_cuda,
            batch_size=batch_size,
            data_and_attributes=data_and_attributes,
            **data_loader_kwargs,
        )
        self.scheme = scheme

    def __iter__(self):
        # TODO: probably don't want to hardcode this condition
        if len(self.labelled_set.indices) == 0 or self.scheme == "alternate":
            return iter(self.full_dataset)
        else:
            return zip(self.full_dataset, cycle(self.labelled_set))
