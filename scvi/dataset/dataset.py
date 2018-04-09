# -*- coding: utf-8 -*-

"""Handling datasets.
For the moment, is initialized with a torch Tensor of size (n_cells, nb_genes)"""
import numpy as np
import torch
from torch.utils.data import Dataset

if torch.cuda.is_available():
    dtype = torch.cuda.FloatTensor
else:
    dtype = torch.FloatTensor


class GeneExpressionDataset(Dataset):
    """Gene Expression dataset. It deals with:
    - log_variational expression -> torch.log(1 + X)
    - local library size normalization (mean, var) per batch
    """

    def __init__(self, Xs):
        # Args:
        # Xs: a list of numpy tensors with .shape[1] identical (total_size*nb_genes)
        self.nb_genes = Xs[0].shape[1]
        self.n_batches = len(Xs)
        assert all(X.shape[1] == self.nb_genes for X in Xs), "All tensors must have same size"

        new_Xs = []
        for i, X in enumerate(Xs):
            X = torch.Tensor(X)
            log_counts = torch.log(torch.sum(X, dim=1))
            local_mean = torch.mean(log_counts)
            local_var = torch.var(log_counts)
            new_Xs += [
                torch.cat((X, local_mean * torch.ones((X.size()[0], 1)),
                           local_var * torch.ones((X.size()[0], 1)),
                           torch.from_numpy(i * np.ones((X.size()[0], 1))).type(torch.FloatTensor)),
                          dim=1)]

        self.X = torch.cat(new_Xs, dim=0)  # .type(dtype)
        self.total_size = self.X.size(0)

    def get_all(self):
        return self.X[:, :-3]

    def get_batches(self):
        dtype = torch.cuda.IntTensor if torch.cuda.is_available() else torch.IntTensor
        return self.X[:, -1].type(dtype).numpy()

    def __len__(self):
        return self.total_size

    def __getitem__(self, idx):
        # Returns the triplet (X, local_mean, local_var)
        return self.X[idx, :self.nb_genes], self.X[idx, self.nb_genes], self.X[idx, self.nb_genes + 1], \
               self.X[idx, self.nb_genes + 2]

    @staticmethod
    def train_test_split(*Xs, train_size=0.75):
        """
        A substitute for the sklearn function to avoid the dependency
        """
        Xs = [np.array(X) for X in Xs]
        split_idx = int(train_size * len(Xs[0]))

        split_list = [[X[:split_idx], X[split_idx:]] for X in Xs]
        return [X for Xs in split_list for X in Xs]
