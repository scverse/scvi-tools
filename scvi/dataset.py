# -*- coding: utf-8 -*-

"""Handling datasets.
For the moment, is initialized with a torch Tensor of size (n_cells, nb_genes)"""
import numpy as np
import torch
from torch.utils.data import Dataset


class GeneExpressionDataset(Dataset):
    """Gene Expression dataset. It deals with:
    - log_variational expression -> torch.log(1 + X)
    - local library size normalization (mean, var) per batch
    """

    def __init__(self, Xs):
        # Args:
        # Xs: a list of numpy tensors with .shape[1] identical (batch_size*nb_genes)
        self.nb_genes = Xs[0].shape[1]
        assert all(X.shape[1] == self.nb_genes for X in Xs), "All tensors must have same size"

        new_Xs = []
        for X in Xs:
            X = torch.Tensor(X)
            log_counts = torch.log(torch.sum(X, dim=1))
            local_mean = torch.mean(log_counts)
            local_var = torch.var(log_counts)
            new_Xs += [
                torch.cat((torch.log(1 + X), local_mean * torch.ones((X.size()[0], 1)),
                           local_var * torch.ones((X.size()[0], 1))),
                          dim=1)]

        self.X = torch.cat(new_Xs, dim=0).type(torch.FloatTensor)
        self.batch_size = self.X.size(0)

    def __len__(self):
        return self.batch_size

    def __getitem__(self, idx):
        # Returns the triplet (X, local_mean, local_var)
        return self.X[idx, :self.nb_genes], self.X[idx, self.nb_genes], self.X[idx, self.nb_genes + 1]


def generate_data(batch_size=20, nb_genes=100):
    # Generating samples according to a ZINB process
    data = np.random.negative_binomial(5, 0.3, size=(batch_size, nb_genes))
    newdata = np.ones((batch_size, nb_genes))
    mask = np.random.binomial(n=1, p=0.7, size=(batch_size, nb_genes))
    for i in range(batch_size):
        newdata[i, :] = data[i, :] / np.sum(data[i, :])
        newdata[i, :] = newdata[i, :] * mask[i, :]
    return newdata
