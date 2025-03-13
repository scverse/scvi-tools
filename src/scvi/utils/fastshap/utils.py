import itertools

import numpy as np
import torch
import torch.nn as nn
from torch.distributions.categorical import Categorical
from torch.utils.data import Dataset


class MaskLayer1d(nn.Module):
    """
    Masking for 1d inputs.

    Args:
      value: replacement value(s) for held out features.
      append: whether to append the mask along feature dimension.
    """

    def __init__(self, value, append):
        super().__init__()
        self.value = value
        self.append = append

    def forward(self, input_tuple):
        """Forward pass."""
        x, S = input_tuple
        x = x * S + self.value * (1 - S)
        if self.append:
            x = torch.cat((x, S), dim=1)
        return x


class KLDivLoss(nn.Module):
    """
    KL divergence loss that applies log softmax operation to predictions.

    Args:
      reduction: how to reduce loss value (e.g., 'batchmean').
      log_target: whether the target is expected as a log probabilities (or as
        probabilities).
    """

    def __init__(self, reduction="batchmean", log_target=False):
        super().__init__()
        self.kld = nn.KLDivLoss(reduction=reduction, log_target=log_target)

    def forward(self, pred, target):
        """
        Evaluate loss.

        Args:
          pred:
          target:
        """
        return self.kld(pred.log_softmax(dim=1), target)


class DatasetRepeat(Dataset):
    """Dataset Repeat.

    Description
    ----------
    A wrapper around multiple datasets that allows repeated elements when the
    dataset sizes don't match. The number of elements is the maximum dataset
    size, and all datasets must be broadcastable to the same size.

    Args:
      datasets: list of dataset objects.

    """

    def __init__(self, datasets):
        # Get maximum number of elements.
        assert np.all([isinstance(dset, Dataset) for dset in datasets])
        items = [len(dset) for dset in datasets]
        num_items = np.max(items)

        # Ensure all datasets align.
        self.dsets = datasets
        self.num_items = num_items
        self.items = items

    def __getitem__(self, index):
        """
        Retrieves the item at the given index.

        Args:
            index (int): The index of the item to retrieve.

        Returns
        -------
            The item at the specified index.
        """
        assert 0 <= index < self.num_items
        return_items = [
            dset[index % num] for dset, num in zip(self.dsets, self.items, strict=False)
        ]
        return tuple(itertools.chain(*return_items))

    def __len__(self):
        """
        Returns the number of items in the collection.

        Returns
        -------
            int: The number of items.
        """
        return self.num_items


class UniformSampler:
    """
    For sampling player subsets with cardinality chosen uniformly at random.

    Args:
      num_players: number of players.
    """

    def __init__(self, num_players):
        self.num_players = num_players

    def sample(self, batch_size):
        """
        Generate sample.

        Args:
          batch_size: number of samples
        """
        rand = torch.rand(batch_size, self.num_players)
        thresh = torch.rand(batch_size, 1)
        S = (thresh > rand).float()

        return S


class ShapleySampler:
    """
    For sampling player subsets from the Shapley distribution.

    Args:
      num_players: number of players.
    """

    def __init__(self, num_players):
        arange = torch.arange(1, num_players)
        w = 1 / (arange * (num_players - arange))
        w = w / torch.sum(w)
        self.categorical = Categorical(probs=w)
        self.num_players = num_players
        self.tril = np.tril(np.ones((num_players - 1, num_players), dtype=np.float32), k=0)
        self.rng = np.random.default_rng()

    def sample(self, batch_size, paired_sampling):
        """
        Generate sample.

        Args:
          batch_size: number of samples.
          paired_sampling: whether to use paired sampling.
        """
        num_included = 1 + self.categorical.sample([batch_size])
        S = self.tril[num_included - 1]
        S = self.rng.permuted(S, axis=1)  # Note: permutes each row.
        if paired_sampling:
            S[1::2] = 1 - S[0 : (batch_size - 1) : 2]  # Note: allows batch_size % 2 == 1.
        return torch.from_numpy(S)
