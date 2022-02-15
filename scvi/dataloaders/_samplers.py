import math
from typing import Iterator, Optional, TypeVar, Union

import numpy as np
import torch
import torch.distributed as dist
from torch.utils.data.distributed import DistributedSampler
from torch.utils.data.sampler import Sampler

T_co = TypeVar("T_co", covariant=True)


class BatchSampler(Sampler):
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
    drop_last
        if int, drops the last batch if its length is less than drop_last.
        if drop_last == True, drops last non-full batch.
        if drop_last == False, iterate over all batches.
    """

    def __init__(
        self,
        indices: np.ndarray,
        batch_size: int,
        shuffle: bool,
        drop_last: Union[bool, int] = False,
    ):
        self.indices = indices
        self.n_obs = len(indices)
        self.batch_size = batch_size
        self.shuffle = shuffle

        if drop_last > batch_size:
            raise ValueError(
                "drop_last can't be greater than batch_size. "
                + "drop_last is {} but batch_size is {}.".format(drop_last, batch_size)
            )

        last_batch_len = self.n_obs % self.batch_size
        if (drop_last is True) or (last_batch_len < drop_last):
            drop_last_n = last_batch_len
        elif (drop_last is False) or (last_batch_len >= drop_last):
            drop_last_n = 0
        else:
            raise ValueError("Invalid input for drop_last param. Must be bool or int.")

        self.drop_last_n = drop_last_n

    def __iter__(self):
        if self.shuffle is True:
            idx = torch.randperm(self.n_obs).tolist()
        else:
            idx = torch.arange(self.n_obs).tolist()

        if self.drop_last_n != 0:
            idx = idx[: -self.drop_last_n]

        data_iter = iter(
            [
                self.indices[idx[i : i + self.batch_size]]
                for i in range(0, len(idx), self.batch_size)
            ]
        )
        return data_iter

    def __len__(self):
        from math import ceil

        if self.drop_last_n != 0:
            length = self.n_obs // self.batch_size
        else:
            length = ceil(self.n_obs / self.batch_size)
        return length


class SubsetDistributedSampler(DistributedSampler):
    """
    Sampler that restricts data loading to a subset of the dataset.

    This wrapper on the PyTorch class adds the ability to choose the subset
    using an iterable of observation-level indexes.

    Parameters
    ----------
    dataset: Dataset used for sampling.
    num_replicas
        Number of processes participating in distributed training. By default,
        :attr:`world_size` is retrieved from the current distributed group.
    rank
        Rank of the current process within :attr:`num_replicas`.
        By default, :attr:`rank` is retrieved from the current distributed
        group.
    shuffle
        If ``True`` (default), sampler will shuffle the indices.
    seed
        random seed used to shuffle the sampler if :attr:`shuffle=True`.
        This number should be identical across all processes in the
        distributed group. Default: ``0``.
    drop_last
        If ``True``, then the sampler will drop the tail of the data to make it
        evenly divisible across the number of replicas. If ``False``, the sampler
        will add extra indices to make the data evenly divisible across the
        replicas. Default: ``False``.
    """

    def __init__(
        self,
        indices: np.ndarray,
        num_replicas: Optional[int] = None,
        rank: Optional[int] = None,
        shuffle: bool = True,
        seed: int = 0,
        drop_last: bool = False,
    ) -> None:

        # PyTorch code only checks length of dataset
        super().__init__(indices, num_replicas, rank, shuffle, seed, drop_last)
        self.indices = np.asarray(indices)

    def __iter__(self) -> Iterator[T_co]:
        iter = super().__iter__(self)
        return iter(torch.from_numpy(self.indices[list(iter)]))
