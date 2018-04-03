# -*- coding: utf-8 -*-

"""Handling datasets.
For the moment, is initialized with a torch Tensor of size (n_cells, nb_genes)"""
import numpy as np
import torch
from torch.utils.data import Dataset


class GeneExpressionDataset(Dataset):
    """Gene Expression dataset."""

    def __init__(self, Xs, dropout=None):
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
                torch.cat((X, local_mean * torch.ones((X.size()[0], 1)), local_var * torch.ones((X.size()[0], 1))),
                          dim=1)]

        self.X = torch.cat(new_Xs, dim=0).type(torch.FloatTensor)
        self.batch_size = self.X.size()[0]
        # self.batch_size = list(X.size())[0]
        # self.nb_genes = list(X.size())[1]
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
            x = self.X[idx, :self.nb_genes] * self.indices_dropout[idx, :]
        else:
            # returns sampled_batch, local_mean, local_var
            x = self.X[idx, :self.nb_genes], self.X[idx, self.nb_genes], self.X[idx, self.nb_genes + 1]
        return x
