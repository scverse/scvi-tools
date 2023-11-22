from torch.utils.data import Dataset, DistributedSampler


class BatchDistributedSampler(DistributedSampler):
    """``EXPERIMENTAL`` Sampler that restricts to loading from a subset of the dataset.

    In contrast to :class:`~torch.utils.data.distributed.DistributedSampler`,
    retrieves a minibatch of data with one call to the dataset's `__getitem__`
    for efficient access to sparse data.

    Parameters
    ----------
    dataset
        :class:`~torch.utils.data.Dataset` instance to sample from.
    batch_size
        Minibatch size to load each iteration for each replica. Thus, the
        effective minibatch size is `batch_size` * `num_replicas`.
    drop_last
        If `True` and the dataset is not evenly divisible by `batch_size`, the last
        incomplete batch is dropped. If `False` and the dataset is not evenly divisible
        by `batch_size`, then the last batch will be smaller than `batch_size`.
    drop_dataset_tail
        If `True` the sampler will drop the tail of the dataset to make it evenly
        divisible by the number of replicas. If `False`, then the sampler will add extra
        indices to make the dataset evenly divisible by the number of replicas.
    **kwargs
        Additional keyword arguments passed into
        :class:`~torch.utils.data.distributed.DistributedSampler`.
    """

    def __init__(
        self,
        dataset: Dataset,
        batch_size: int = 128,
        drop_last: bool = False,
        drop_dataset_tail: bool = False,
        **kwargs,
    ):
        super().__init__(dataset, drop_last=drop_dataset_tail, **kwargs)
        self.batch_size = batch_size
        self.drop_last_batch = drop_last  # drop_last already defined in parent

    def __iter__(self):
        """Iterates over a subset of indices from the dataset.

        Based on :meth:`~torch.utils.data.BatchSampler.__iter__`.

        Notes
        -----
        `super().__iter__()` iterates over a subset of the dataset based on the current
        `rank` and `num_replicas`.
        """
        sampler_iter = super().__iter__()
        if self.drop_last_batch:
            while True:
                try:
                    batch = [next(sampler_iter) for _ in range(self.batch_size)]
                    yield batch
                except StopIteration:
                    break
        else:
            batch = [0] * self.batch_size
            idx_in_batch = 0
            for idx in sampler_iter:
                batch[idx_in_batch] = idx
                idx_in_batch += 1
                if idx_in_batch == self.batch_size:
                    yield batch
                    idx_in_batch = 0
                    batch = [0] * self.batch_size
            if idx_in_batch > 0:
                yield batch[:idx_in_batch]
