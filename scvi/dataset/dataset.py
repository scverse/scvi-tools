# -*- coding: utf-8 -*-

"""Handling datasets.
For the moment, is initialized with a torch Tensor of size (n_cells, nb_genes)"""
import os

import numpy as np
import scipy.sparse as sp_sparse
import torch
import urllib.request
from torch.utils.data import Dataset


class GeneExpressionDataset(Dataset):
    """Gene Expression dataset. It deals with:
    - log_variational expression -> torch.log(1 + X)
    - local library size normalization (mean, var) per batch
    """

    def __init__(self, X, local_means, local_vars, batch_indices, labels, gene_names=None, n_batches=1):
        # Args:
        # Xs: a list of numpy tensors with .shape[1] identical (total_size*nb_genes)
        # or a list of scipy CSR sparse matrix,
        # or transposed CSC sparse matrix (the argument sparse must then be set to true)
        self.total_size, self.nb_genes = X.shape
        self.n_batches = n_batches
        self.local_means = local_means
        self.local_vars = local_vars
        self.batch_indices = batch_indices
        self.X = X
        self.labels = labels

        if gene_names is not None:
            self.gene_names = np.char.upper(gene_names)  # Take an upper case convention for gene names

    def __len__(self):
        return self.total_size

    def __getitem__(self, idx):
        return self.X[idx].toarray()[0], self.local_means[idx], self.local_vars[idx], \
               self.batch_indices[idx], self.labels[idx]

    def download(self):
        r = urllib.request.urlopen(self.url)
        print("Downloading data")

        def readIter(f, blocksize=1000):
            """Given a file 'f', returns an iterator that returns bytes of
            size 'blocksize' from the file, using read()."""
            while True:
                data = f.read(blocksize)
                if not data:
                    break
                yield data

        # Create the path to save the data
        if not os.path.exists(self.save_path):
            os.makedirs(self.save_path)

        with open(self.save_path + self.download_name, 'wb') as f:
            for data in readIter(r):
                f.write(data)

    def download_and_preprocess(self):
        if not os.path.exists(self.save_path + self.download_name):
            self.download()
        return self.preprocess()

    @staticmethod
    def train_test_split(*Xs, train_size=0.75):
        """
        A substitute for the sklearn function to avoid the dependency
        """
        Xs = [np.array(X) for X in Xs]
        split_idx = int(train_size * len(Xs[0]))
        all_indices = np.arange(len(Xs[0]))
        train_indices = np.random.choice(all_indices, size=split_idx, replace=False)
        test_indices = np.array(list(set(all_indices).difference(set(train_indices))))
        split_list = [[X[train_indices], X[test_indices]] for X in Xs]
        return [X for Xs in split_list for X in Xs]

    @staticmethod
    def get_attributes_from_matrix(X, batch_index=0, labels=None):
        log_counts = torch.log(torch.from_numpy(np.asarray((X.sum(axis=1)), dtype=float)).type(torch.FloatTensor))
        local_mean = torch.mean(log_counts) * torch.ones((X.shape[0], 1))
        local_var = torch.var(log_counts) * torch.ones((X.shape[0], 1))
        batch_index = torch.LongTensor(batch_index * np.ones((X.shape[0], 1)))
        labels = torch.LongTensor(labels.reshape(-1, 1)) if labels is not None else torch.zeros_like(batch_index)
        return X, local_mean, local_var, batch_index, labels

    @staticmethod
    def get_attributes_from_list(Xs, list_labels=None):
        nb_genes = Xs[0].shape[1]
        # n_batches = len(Xs)
        assert all(X.shape[1] == nb_genes for X in Xs), "All tensors must have same size"

        new_Xs = []
        local_means = []
        local_vars = []
        batch_indices = []
        labels = []
        for i, X in enumerate(Xs):
            label = list_labels[i] if list_labels is not None else list_labels
            X, local_mean, local_var, batch_index, label = (
                GeneExpressionDataset.get_attributes_from_matrix(X, batch_index=i, labels=label)
            )
            new_Xs += [X]
            local_means += [local_mean]
            local_vars += [local_var]
            batch_indices += [batch_index]
            labels += [label]
        local_means = torch.cat(local_means)
        local_vars = torch.cat(local_vars)
        batch_indices = torch.cat(batch_indices)
        labels = torch.cat(labels)
        X = sp_sparse.vstack(new_Xs)
        return X, local_means, local_vars, batch_indices, labels
