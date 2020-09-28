import copy
import logging
from typing import Dict, List, Optional, Union

import anndata
import numpy as np
import torch
from torch.utils.data import DataLoader

from scvi import _CONSTANTS
from scvi.core._log_likelihood import (
    compute_elbo,
    compute_marginal_log_likelihood_autozi,
    compute_marginal_log_likelihood_scvi,
    compute_reconstruction_error,
)
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

    def __init__(self, indices: np.ndarray, batch_size: int, shuffle: bool):
        self.indices = indices
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
        return len(self.indices) // self.batch_size


class ScviDataLoader:
    """
    Scvi Data Loader.

    A `ScviDataLoader` instance is instantiated with a model and a gene_dataset, and
    as well as additional arguments that for Pytorch's `DataLoader`. A subset of indices can be specified, for
    purposes such as splitting the data into train/test or labelled/unlabelled (for semi-supervised learning).
    Each trainer instance of the `Trainer` class can therefore have multiple `ScviDataLoader` instances to train a model.
    A `ScviDataLoader` instance also comes with methods to compute likelihood and other relevant training metrics.

    Parameters
    ----------
    model
        A model instance from class ``VAE``, ``VAEC``, ``SCANVI``
    gene_dataset
        A gene_dataset instance like ``CortexDataset()``
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
        model,
        adata: anndata.AnnData,
        shuffle=False,
        indices=None,
        use_cuda=True,
        batch_size=128,
        data_loader_kwargs=dict(),
    ):
        self.model = model
        if "_scvi" not in adata.uns.keys():
            raise ValueError("Please run setup_anndata() on your anndata object first.")
        for key in self._data_and_attributes.keys():
            if key not in adata.uns["_scvi"]["data_registry"].keys():
                raise ValueError(
                    "{} required for model but not included when setup_anndata was run".format(
                        key
                    )
                )
        self.dataset = ScviDataset(adata, getitem_tensors=self._data_and_attributes)
        self.to_monitor = []
        self.use_cuda = use_cuda

        if indices is None:
            inds = np.arange(len(self.dataset))
            if shuffle:
                sampler_kwargs = {
                    "indices": inds,
                    "batch_size": batch_size,
                    "shuffle": True,
                }
            else:
                sampler_kwargs = {
                    "indices": inds,
                    "batch_size": batch_size,
                    "shuffle": False,
                }
        else:
            if hasattr(indices, "dtype") and indices.dtype is np.dtype("bool"):
                indices = np.where(indices)[0].ravel()
            indices = np.asarray(indices)
            sampler_kwargs = {
                "indices": indices,
                "batch_size": batch_size,
                "shuffle": True,
            }

        self.sampler_kwargs = sampler_kwargs
        sampler = BatchSampler(**self.sampler_kwargs)
        self.data_loader_kwargs = copy.copy(data_loader_kwargs)
        # do not touch batch size here, sampler gives batched indices
        self.data_loader_kwargs.update({"sampler": sampler, "batch_size": None})
        self.data_loader = DataLoader(self.dataset, **self.data_loader_kwargs)
        self.original_indices = self.indices

    @property
    def _data_and_attributes(self):
        """
        Returns dictionary where key is scVI data field and value is its datatype.

        Used by ScviDataset when _get_item_ is called by mapping the keys in the dict
        to its associated datatype
        """
        return {
            _CONSTANTS.X_KEY: np.float32,
            _CONSTANTS.BATCH_KEY: np.int64,
            _CONSTANTS.LOCAL_L_MEAN_KEY: np.float32,
            _CONSTANTS.LOCAL_L_VAR_KEY: np.float32,
            _CONSTANTS.LABELS_KEY: np.int64,
        }

    def accuracy(self):
        pass

    accuracy.mode = "max"

    @property
    def indices(self) -> np.ndarray:
        """Returns the current dataloader indices used by the object."""
        if hasattr(self.data_loader.sampler, "indices"):
            return self.data_loader.sampler.indices
        else:
            return np.arange(len(self.dataset))

    @property
    def n_cells(self) -> int:
        """Returns the number of studied cells."""
        if hasattr(self.data_loader.sampler, "indices"):
            return len(self.data_loader.sampler.indices)
        else:
            return self.dataset.n_cells

    @property
    def scvi_data_loader_type(self) -> str:
        """Returns the dataloader class name."""
        return self.__class__.__name__

    def __iter__(self):
        return map(self.to_cuda, iter(self.data_loader))

    def to_cuda(self, tensors: Dict[str, torch.Tensor]) -> Dict[str, torch.Tensor]:
        """
        Converts dict of tensors to cuda.

        Parameters
        ----------
        tensors
            tensors to convert

        """
        return {k: (t.cuda() if self.use_cuda else t) for k, t in tensors.items()}

    def update(self, data_loader_kwargs: dict) -> "ScviDataLoader":
        """
        Updates the dataloader.

        Parameters
        ----------
        data_loader_kwargs
            dataloader updates.

        Returns
        -------
        Updated ScviDataLoader

        """
        scdl = copy.copy(self)
        scdl.data_loader_kwargs = copy.copy(self.data_loader_kwargs)
        scdl.data_loader_kwargs.update(data_loader_kwargs)
        scdl.data_loader = DataLoader(self.dataset, **scdl.data_loader_kwargs)
        return scdl

    def update_batch_size(self, batch_size):
        self.sampler_kwargs.update({"batch_size": batch_size})
        sampler = BatchSampler(**self.sampler_kwargs)
        return self.update({"sampler": sampler, "batch_size": None})

    def sequential(self, batch_size: Optional[int] = 128) -> "ScviDataLoader":
        """
        Returns a copy of the object that iterate over the data sequentially.

        Parameters
        ----------
        batch_size
            New batch size.

        """
        self.sampler_kwargs = {
            "indices": self.indices,
            "batch_size": batch_size,
            "shuffle": False,
        }
        return self.update({"sampler": BatchSampler(**self.sampler_kwargs)})

    @torch.no_grad()
    def elbo(self) -> torch.Tensor:
        """Returns the Evidence Lower Bound associated to the object."""
        elbo = compute_elbo(self.model, self)
        logger.debug("ELBO : %.4f" % elbo)
        return elbo

    elbo.mode = "min"

    @torch.no_grad()
    def reconstruction_error(self) -> torch.Tensor:
        """Returns the reconstruction error associated to the object."""
        reconstruction_error = compute_reconstruction_error(self.model, self)
        logger.debug("Reconstruction Error : %.4f" % reconstruction_error)
        return reconstruction_error

    reconstruction_error.mode = "min"

    @torch.no_grad()
    def marginal_ll(self, n_mc_samples: Optional[int] = 1000) -> torch.Tensor:
        """
        Estimates the marginal likelihood of the object's data.

        Parameters
        ----------
        n_mc_samples
            Number of MC estimates to use

        Returns
        -------
        Marginal LL

        """
        if (
            hasattr(self.model, "reconstruction_loss")
            and self.model.reconstruction_loss == "autozinb"
        ):
            ll = compute_marginal_log_likelihood_autozi(self.model, self, n_mc_samples)
        else:
            ll = compute_marginal_log_likelihood_scvi(self.model, self, n_mc_samples)
        logger.debug("True LL : %.4f" % ll)
        return ll

    def update_sampler_indices(self, idx: Union[List, np.ndarray]):
        """
        Updates the data loader indices.

        More precisely, this method can be used to temporarily change which cells __iter__
        will yield. This is particularly useful for computational considerations when one is only interested
        in a subset of the cells of the ScviDataLoader object.
        This method should be used carefully and requires to reset the data loader to its
        original value after use.

        Parameters
        ----------
        idx :
            Indices (in [0, len(dataset)] to sample from

        Examples
        --------
        >>> old_loader = self.data_loader
        >>> cell_indices = np.array([1, 2, 3])
        >>> self.update_sampler_indices(cell_indices)
        >>> for tensors in self:
        >>>    # your code

        >>> # Do not forget next line!
        >>> self.data_loader = old_loader

        """
        self.sampler_kwargs.update({"indices": idx})
        sampler = BatchSampler(**self.sampler_kwargs)
        self.data_loader_kwargs.update({"sampler": sampler, "batch_size": None})
        self.data_loader = DataLoader(self.dataset, **self.data_loader_kwargs)
