# -*- coding: utf-8 -*-

"""Handling datasets.
For the moment, is initialized with a torch Tensor of size (n_cells, nb_genes)"""
import os

import numpy as np
import scipy.sparse as sp_sparse
import torch
import urllib.request
from sklearn.preprocessing import StandardScaler
from torch.utils.data import Dataset


class GeneExpressionDataset(Dataset):
    """Gene Expression dataset. It deals with:
    - log_variational expression -> torch.log(1 + X)
    - local library size normalization (mean, var) per batch
    """

    def __init__(self, X, local_means, local_vars, batch_indices, labels, gene_names=None):
        # Args:
        # Xs: a list of numpy tensors with .shape[1] identical (total_size*nb_genes)
        # or a list of scipy CSR sparse matrix,
        # or transposed CSC sparse matrix (the argument sparse must then be set to true)
        self.X = X
        self.total_size, self.nb_genes = self.X.shape
        self.local_means = local_means
        self.local_vars = local_vars
        self.labels = labels
        self.batch_indices = batch_indices
        self.n_batches = len(np.unique(batch_indices))
        self.dense = type(self.X) is np.ndarray
        self.n_labels = len(np.unique(labels))

        if gene_names is not None:
            self.gene_names = np.char.upper(gene_names)  # Take an upper case convention for gene names

    def __len__(self):
        return self.X.shape[0]

    def __getitem__(self, idx):
        return idx

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

    def collate_fn(self, batch):
        indexes = np.array(batch)
        X = torch.FloatTensor(self.X[indexes]) if self.dense else torch.FloatTensor(self.X[indexes].toarray())
        return X, torch.FloatTensor(self.local_means[indexes]), \
            torch.FloatTensor(self.local_vars[indexes]), \
            torch.LongTensor(self.batch_indices[indexes]), \
            torch.LongTensor(self.labels[indexes])

    def subsample_genes(self, p_genes=1., subset_genes=None):
        n_cells, n_genes = self.X.shape
        new_n_genes = int(p_genes * n_genes) if type(p_genes) is not int else p_genes
        if subset_genes is None:
            print("Downsampling from %i to %i genes" % (n_genes, new_n_genes))
            std_scaler = StandardScaler(with_mean=False)
            std_scaler.fit(self.X.astype(np.float64))
            subset_genes = np.argsort(std_scaler.var_)[::-1][:new_n_genes]
        else:
            print("Downsampling from %i to %i genes" % (n_genes, len(subset_genes)))
        self.X = self.X[:, subset_genes]
        self.nb_genes = self.X.shape[1]
        self.gene_names = self.gene_names[subset_genes]

    def subsample_cells(self, p_cells=1.):
        n_cells, n_genes = self.X.shape
        new_n_cells = int(p_cells * n_genes) if type(p_cells) is not int else p_cells
        print("Downsampling from %i to %i cells" % (n_cells, new_n_cells))
        indices = np.argsort(self.X.sum(axis=1))[::-1][:new_n_cells]

        def _subsample(a, indices):
            return a[indices]

        self.X, self.local_means, self.local_vars, self.batch_indices, self.labels = (
            _subsample(a, indices) for a in (self.X, self.local_means, self.local_vars,
                                             self.batch_indices, self.labels))

    @staticmethod
    def get_attributes_from_matrix(X, batch_index=0, labels=None):
        log_counts = np.log(X.sum(axis=1))
        local_mean = (np.mean(log_counts) * np.ones((X.shape[0], 1))).astype(np.float32)
        local_var = (np.var(log_counts) * np.ones((X.shape[0], 1))).astype(np.float32)
        batch_index = batch_index * np.ones((X.shape[0], 1))
        labels = labels.reshape(-1, 1) if labels is not None else np.zeros_like(batch_index)
        return X, local_mean, local_var, batch_index, labels

    @staticmethod
    def get_attributes_from_list(Xs, list_labels=None):
        nb_genes = Xs[0].shape[1]
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
        local_means = np.concatenate(local_means)
        local_vars = np.concatenate(local_vars)
        batch_indices = np.concatenate(batch_indices)
        labels = np.concatenate(labels)
        X = np.concatenate(new_Xs) if type(new_Xs[0]) is np.ndarray else sp_sparse.vstack(new_Xs)
        return X, local_means, local_vars, batch_indices, labels
