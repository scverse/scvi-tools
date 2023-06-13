from torch.utils.data import DistributedSampler


class BatchDistributedSampler(DistributedSampler):
    """Sampler that restricts data loading to a subset of the dataset.

    In contrast to :class:`~torch.utils.data.distributed.DistributedSampler`, this
    sampler retrieves a minibatch of data with one call to the dataset's
    `__getitem__` for efficient access to sparse matrices.
    """

    def __init__(self, *args, batch_size: int = 128, **kwargs):
        super().__init__(*args, **kwargs)
        self.batch_size = batch_size

    def __iter__(self):
        """Based on :meth:`~torch.utils.data.BatchSampler.__iter__`."""
        if self.drop_last:
            sampler_iter = super().__iter__()
            while True:
                try:
                    batch = [next(sampler_iter) for _ in range(self.batch_size)]
                    yield batch
                except StopIteration:
                    break
        else:
            batch = [0] * self.batch_size
            idx_in_batch = 0
            for idx in super().__iter__():
                batch[idx_in_batch] = idx
                idx_in_batch += 1
                if idx_in_batch == self.batch_size:
                    yield batch
                    idx_in_batch = 0
                    batch = [0] * self.batch_size
            if idx_in_batch > 0:
                yield batch[:idx_in_batch]
