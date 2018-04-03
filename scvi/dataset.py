# -*- coding: utf-8 -*-

"""Handling datasets.
For the moment, is initialized with a torch Tensor of size (n_cells, nb_genes)"""
from torch.utils.data import Dataset


class GeneExpressionDataset(Dataset):
    """Gene Expression dataset."""

    def __init__(self, X):

        # Args:
        # X: torch tensor of size (batch_size*nb_genes)

        self.X = X
        self.batch_size = list(X.size())[0]
        self.nb_genes = list(X.size())[1]

    def __len__(self):
        return self.batch_size

    def __getitem__(self, idx):
        x = self.X[idx, :]
        return x
