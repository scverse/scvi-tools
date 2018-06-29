# -*- coding: utf-8 -*-

"""Handling datasets.
For the moment, is initialized with a torch Tensor of size (n_cells, nb_genes)"""
import os
import urllib.request

import numpy as np
import scipy.sparse as sp_sparse
import torch
from sklearn.preprocessing import StandardScaler
from torch.utils.data import Dataset

from scvi.dataset.utils import filter_genes, arrange_categories


class GeneExpressionDataset(Dataset):
    """Gene Expression dataset. It deals with:
    - log_variational expression -> torch.log(1 + X)
    - local library size normalization (mean, var) per batch
    """

    def __init__(self, X, local_means, local_vars, batch_indices, labels, gene_names=None, cell_types=None):
        # Args:
        # Xs: a list of numpy tensors with .shape[1] identical (total_size*nb_genes)
        # or a list of scipy CSR sparse matrix,
        # or transposed CSC sparse matrix (the argument sparse must then be set to true)
        self.X = X
        self.nb_genes = self.X.shape[1]
        self.dense = type(self.X) is np.ndarray
        self.local_means = local_means
        self.local_vars = local_vars
        self.batch_indices, self.n_batches = arrange_categories(batch_indices)
        self.labels, self.n_labels = arrange_categories(labels)

        if gene_names is not None:
            assert self.nb_genes == len(gene_names)
            self.gene_names = np.array(gene_names, dtype=np.str)
        if cell_types is not None:
            assert self.n_labels == len(cell_types)
            self.cell_types = np.array(cell_types, dtype=np.str)

    def __len__(self):
        return self.X.shape[0]

    def __getitem__(self, idx):
        return idx

    def download_and_preprocess(self):
        self.download()
        return self.preprocess()

    def collate_fn(self, batch):
        indexes = np.array(batch)
        X = torch.FloatTensor(self.X[indexes]) if self.dense else torch.FloatTensor(self.X[indexes].toarray())
        return X, torch.FloatTensor(self.local_means[indexes]), \
            torch.FloatTensor(self.local_vars[indexes]), \
            torch.LongTensor(self.batch_indices[indexes]), \
            torch.LongTensor(self.labels[indexes])

    def update_genes(self, subset_genes):
        if hasattr(self, 'gene_names'):
            self.gene_names = self.gene_names[subset_genes]
        if hasattr(self, 'gene_symbols'):
            self.gene_symbols = self.gene_symbols[subset_genes]
        self.nb_genes = self.X.shape[1]

    def update_cells(self, subset_cells):
        new_n_cells = len(subset_cells) if subset_cells.dtype is not np.bool else subset_cells.sum()
        print("Downsampling from %i to %i cells" % (len(self), new_n_cells))
        for attr_name in ['X', 'local_means', 'local_vars', 'labels', 'batch_indices']:
            setattr(self, attr_name, getattr(self, attr_name)[subset_cells])

    def subsample_genes(self, new_n_genes=None, subset_genes=None):
        n_cells, n_genes = self.X.shape
        if subset_genes is None and \
                (not hasattr(self, 'gene_names') or new_n_genes is False or new_n_genes >= n_genes):
            return None  # Do nothing if subsample more genes than total number of genes
        if subset_genes is None:
            print("Downsampling from %i to %i genes" % (n_genes, new_n_genes))
            std_scaler = StandardScaler(with_mean=False)
            std_scaler.fit(self.X.astype(np.float64))
            subset_genes = np.argsort(std_scaler.var_)[::-1][:new_n_genes]
        else:
            print("Downsampling from %i to %i genes" % (n_genes, len(subset_genes)))
        self.X = self.X[:, subset_genes]
        self.update_genes(subset_genes)

    def filter_genes(self, gene_names_ref, on='gene_names'):
        """
        Same as filter_genes but overwrites on current dataset instead of returning data,
        and updates genes names and symbols
        """
        self.X, subset_genes = filter_genes(self, gene_names_ref, on=on)
        self.update_genes(subset_genes)

    def subsample_cells(self, size=1.):
        n_cells, n_genes = self.X.shape
        new_n_cells = int(size * n_genes) if type(size) is not int else size
        indices = np.argsort(self.X.sum(axis=1))[::-1][:new_n_cells]
        self.update_cells(indices)

    def filter_cell_types(self, cell_types):
        """
        :param cell_types: numpy array of type np.int (indices) or np.str (cell-types names)
        :return:
        """
        if type(cell_types[0]) is not int:
            current_cell_types = list(self.cell_types)
            cell_types_idx = np.array([current_cell_types.index(cell_type) for cell_type in cell_types])
        else:
            cell_types_idx = np.array(cell_types, dtype=np.int64)
        if hasattr(self, 'cell_types'):
            self.cell_types = self.cell_types[cell_types_idx]
            print("Only keeping cell types: \n" + '\n'.join(list(self.cell_types)))
        idx_to_keep = []
        for idx in cell_types_idx:
            idx_to_keep += [np.where(self.labels == idx)[0]]
        self.update_cells(np.concatenate(idx_to_keep))
        self.labels, self.n_labels = arrange_categories(self.labels, mapping_from=cell_types_idx)

    def download(self):
        if hasattr(self, 'urls') and hasattr(self, 'download_names'):
            for url, download_name in zip(self.urls, self.download_names):
                GeneExpressionDataset._download(url, self.save_path, download_name)
        elif hasattr(self, 'url') and hasattr(self, 'download_name'):
            GeneExpressionDataset._download(self.url, self.save_path, self.download_name)

    @staticmethod
    def _download(url, save_path, download_name):
        if os.path.exists(save_path + download_name):
            print("File %s already downloaded" % (save_path + download_name))
            return
        r = urllib.request.urlopen(url)
        print("Downloading file at %s" % save_path + download_name)

        def readIter(f, blocksize=1000):
            """Given a file 'f', returns an iterator that returns bytes of
            size 'blocksize' from the file, using read()."""
            while True:
                data = f.read(blocksize)
                if not data:
                    break
                yield data

        # Create the path to save the data
        if not os.path.exists(save_path):
            os.makedirs(save_path)

        with open(save_path + download_name, 'wb') as f:
            for data in readIter(r):
                f.write(data)

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

    @staticmethod
    def concat_datasets(*gene_datasets, on='gene_names', shared_labels=True):
        """
        Combines multiple unlabelled gene_datasets based on the intersection of gene names intersection.
        Datasets should all have gene_dataset.n_labels=0.
        Batch indices are generated in the same order as datasets are given.
        :param gene_datasets: a sequence of gene_datasets object
        :return: a GeneExpressionDataset instance of the concatenated datasets
        """
        assert all([hasattr(gene_dataset, on) for gene_dataset in gene_datasets])

        gene_names_ref = set.intersection(*[set(getattr(gene_dataset, on)) for gene_dataset in gene_datasets])
        # keep gene order of the first dataset
        gene_names_ref = [gene_name for gene_name in getattr(gene_datasets[0], on) if gene_name in gene_names_ref]
        print("Keeping %d genes" % len(gene_names_ref))

        Xs = [filter_genes(gene_dataset, gene_names_ref, on=on)[0] for gene_dataset in gene_datasets]
        if gene_datasets[0].dense:
            X = np.concatenate([X if type(X) is np.ndarray else X.A for X in Xs])
        else:
            X = sp_sparse.vstack([X if type(X) is not np.ndarray else sp_sparse.csr_matrix(X) for X in Xs])

        batch_indices = np.zeros((X.shape[0], 1))
        n_batch_offset = 0
        current_index = 0
        for gene_dataset in gene_datasets:
            next_index = current_index + len(gene_dataset)
            batch_indices[current_index:next_index] = gene_dataset.batch_indices + n_batch_offset
            n_batch_offset += gene_dataset.n_batches
            current_index = next_index

        cell_types = None
        if shared_labels:
            if all([hasattr(gene_dataset, "cell_types") for gene_dataset in gene_datasets]):
                cell_types = list(
                    set([cell_type for gene_dataset in gene_datasets for cell_type in gene_dataset.cell_types])
                )
                labels = []
                for gene_dataset in gene_datasets:
                    mapping = [cell_types.index(cell_type) for cell_type in gene_dataset.cell_types]
                    labels += [arrange_categories(gene_dataset.labels, mapping_to=mapping)[0]]
                labels = np.concatenate(labels)
            else:
                labels = np.concatenate([gene_dataset.labels for gene_dataset in gene_datasets])
        else:
            labels = np.zeros((X.shape[0], 1))
            n_labels_offset = 0
            current_index = 0
            for gene_dataset in gene_datasets:
                next_index = current_index + len(gene_dataset)
                labels[current_index:next_index] = gene_dataset.labels + n_labels_offset
                n_labels_offset += gene_dataset.n_labels
                current_index = next_index

        _, local_mean, local_var, _, _ = GeneExpressionDataset.get_attributes_from_matrix(X)
        return GeneExpressionDataset(X, local_mean, local_var, batch_indices, labels,
                                     gene_names=gene_names_ref, cell_types=cell_types)
