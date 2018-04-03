# -*- coding: utf-8 -*-

"""Handling datasets.
For the moment, is initialized with a torch Tensor of size (n_cells, nb_genes)"""
import numpy as np
import torch
from torch.utils.data import Dataset


class GeneExpressionDataset(Dataset):
    """Gene Expression dataset."""

    def __init__(self, X, dropout=None):

        # Args:
        # X: torch tensor of size (batch_size*nb_genes)

        self.X = X
        self.batch_size = list(X.size())[0]
        self.nb_genes = list(X.size())[1]
        self.dropout = dropout
        if self.dropout is not None:
            self.indices_dropout = torch.FloatTensor(
                np.random.binomial(1, self.dropout, size=(self.batch_size, self.nb_genes)))

    def __len__(self):
        return self.batch_size

    def __getitem__(self, idx):
        if self.dropout is not None:
            # 1. Generate the dropout indices on the fly
            # x = self.X[idx, :]*torch.FloatTensor(np.random.binomial(1,self.dropout,size=self.nb_genes))

            # 2. Use previously generated indices
            x = self.X[idx, :] * self.indices_dropout[idx, :]
        else:
            x = self.X[idx, :]
        return x
