import copy
import logging

import numpy as np
from torch.utils.data import (
    BatchSampler,
    DataLoader,
    RandomSampler,
    Sampler,
    SequentialSampler,
)

from scvi import settings
from scvi.data import AnnDataManager

from ._samplers import BatchDistributedSampler
from ._semi_dataloader import labelled_indices_generator, subsample_labels

logger = logging.getLogger(__name__)


class AnnDataLoader(DataLoader):
    """DataLoader for loading tensors from AnnData objects.

    Parameters
    ----------
    adata_manager
        :class:`~scvi.data.AnnDataManager` object with a registered AnnData object.
    indices
        The indices of the observations in `adata_manager.adata` to load.
    batch_size
        Minibatch size to load each iteration. If `distributed_sampler` is `True`,
        refers to the minibatch size per replica. Thus, the effective minibatch
        size is `batch_size` * `num_replicas`.
    shuffle
        Whether the dataset should be shuffled prior to sampling.
    sampler
        Defines the strategy to draw samples from the dataset. Can be any Iterable with __len__
        implemented. If specified, shuffle must not be specified. By default, we use a custom
        sampler that is designed to get a minibatch of data with one call to __getitem__.
    drop_last
        If `True` and the dataset is not evenly divisible by `batch_size`, the last
        incomplete batch is dropped. If `False` and the dataset is not evenly divisible
        by `batch_size`, then the last batch will be smaller than `batch_size`.
    drop_dataset_tail
        Only used if `distributed_sampler` is `True`. If `True` the sampler will drop
        the tail of the dataset to make it evenly divisible by the number of replicas.
        If `False`, then the sampler will add extra indices to make the dataset evenly
        divisible by the number of replicas.
    data_and_attributes
        Dictionary with keys representing keys in data registry (``adata_manager.data_registry``)
        and value equal to desired numpy loading type (later made into torch tensor) or list of
        such keys. A list can be used to subset to certain keys in the event that more tensors than
        needed have been registered. If ``None``, defaults to all registered data.
    iter_ndarray
        Whether to iterate over numpy arrays instead of torch tensors
    distributed_sampler
        ``EXPERIMENTAL`` Whether to use :class:`~scvi.dataloaders.BatchDistributedSampler` as the
        sampler. If `True`, `sampler` must be `None`.
    n_samples_per_label
        Number of subsamples for each label class to sample per epoch
    load_sparse_tensor
        ``EXPERIMENTAL`` If ``True``, loads data with sparse CSR or CSC layout as a
        :class:`~torch.Tensor` with the same layout. Can lead to speedups in data transfers to
        GPUs, depending on the sparsity of the data.
    **kwargs
        Additional keyword arguments passed into :class:`~torch.utils.data.DataLoader`.

    Notes
    -----
    If `sampler` is not specified, a :class:`~torch.utils.data.BatchSampler` instance is
    passed in as the sampler, which retrieves a minibatch of data with one call to
    :meth:`~scvi.data.AnnTorchDataset.__getitem__`. This is useful for fast access to
    sparse matrices as retrieving single observations and then collating is inefficient.
    """

    def __init__(
        self,
        adata_manager: AnnDataManager,
        indices: list[int] | list[bool] | None = None,
        batch_size: int = 128,
        shuffle: bool = False,
        sampler: Sampler | None = None,
        drop_last: bool = False,
        drop_dataset_tail: bool = False,
        data_and_attributes: list[str] | dict[str, np.dtype] | None = None,
        iter_ndarray: bool = False,
        distributed_sampler: bool = False,
        n_samples_per_label: int | None = None,
        load_sparse_tensor: bool = False,
        **kwargs,
    ):
        if indices is None:
            indices = np.arange(adata_manager.adata.shape[0])
        else:
            if hasattr(indices, "dtype") and indices.dtype is np.dtype("bool"):
                indices = np.where(indices)[0].ravel()
            indices = np.asarray(indices)
        self.indices = indices
        self.n_samples_per_label = n_samples_per_label
        self.dataset = adata_manager.create_torch_dataset(
            indices=indices,
            data_and_attributes=data_and_attributes,
            load_sparse_tensor=load_sparse_tensor,
        )
        if "num_workers" not in kwargs:
            kwargs["num_workers"] = settings.dl_num_workers
        if "persistent_workers" not in kwargs:
            kwargs["persistent_workers"] = settings.dl_persistent_workers

        self.kwargs = copy.deepcopy(kwargs)
        self.adata_manager = adata_manager
        self.data_and_attributes = data_and_attributes
        self._shuffle = shuffle
        self._batch_size = batch_size
        self._drop_last = drop_last
        self.load_sparse_tensor = load_sparse_tensor

        if sampler is not None and distributed_sampler:
            raise ValueError("Cannot specify both `sampler` and `distributed_sampler`.")

        # Next block of code is for the case of labeled anndataloder used in scanvi multigpu:
        self.labeled_locs, labelled_idx = [], []
        if adata_manager.registry["model_name"] == "SCANVI":
            self.labeled_locs, labelled_idx = labelled_indices_generator(
                adata_manager, indices, self.indices, self.n_samples_per_label
            )

        # custom sampler for efficient minibatching on sparse matrices
        if sampler is None:
            if not distributed_sampler:
                sampler_cls = SequentialSampler if not shuffle else RandomSampler
                sampler = BatchSampler(
                    sampler=sampler_cls(self.dataset),
                    batch_size=batch_size,
                    drop_last=drop_last,
                )
            else:
                sampler = BatchDistributedSampler(
                    self.dataset,
                    batch_size=batch_size,
                    drop_last=drop_last,
                    drop_dataset_tail=drop_dataset_tail,
                    shuffle=shuffle,
                )
            # do not touch batch size here, sampler gives batched indices
            # This disables PyTorch automatic batching, which is necessary
            # for fast access to sparse matrices
            self.kwargs.update({"batch_size": None, "shuffle": False})

        self.kwargs.update({"sampler": sampler})

        if iter_ndarray:
            self.kwargs.update({"collate_fn": lambda x: x})

        super().__init__(self.dataset, **self.kwargs)

    def resample_labels(self):
        """Resamples the labeled data."""
        self.kwargs.pop("batch_size", None)
        self.kwargs.pop("shuffle", None)
        self.kwargs.pop("sampler", None)
        self.kwargs.pop("collate_fn", None)
        AnnDataLoader(
            self.adata_manager,
            indices=subsample_labels(self.labeled_locs, self.n_samples_per_label),
            shuffle=self._shuffle,
            batch_size=self._batch_size,
            data_and_attributes=self.data_and_attributes,
            drop_last=self._drop_last,
            **self.kwargs,
        )
