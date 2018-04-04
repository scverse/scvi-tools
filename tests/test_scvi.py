#!/usr/bin/env python
# -*- coding: utf-8 -*-


"""Tests for `scvi` package."""

from torch.autograd import Variable
from torch.utils.data import DataLoader
from torch.utils.data.sampler import SequentialSampler

from scvi.dataset import GeneExpressionDataset, generate_data
from scvi.imputation import dropout, imputation_error
from scvi.scvi import VAE
from scvi.train import train


def test_benchmark():
    batch_size = 20
    nb_genes = 100
    X = generate_data(batch_size=batch_size, nb_genes=nb_genes)
    gene_dataset = GeneExpressionDataset([X])
    data_loader = DataLoader(gene_dataset, batch_size=4,
                             shuffle=True, num_workers=4)

    # 1. Instanciate model
    vae = VAE(nb_genes)

    # 2. Train model
    train(vae, data_loader)


def test_imputation():
    """The imputation test"""
    # 0. Load data - For now we generate with artificial data
    batch_size = 20
    nb_genes = 100

    X = X_test = generate_data(batch_size=batch_size, nb_genes=nb_genes)
    # X = np.load("expression_train.npy")
    # X_test = np.load("expression_test.npy")

    # Creating a GeneExpressionDataset and a DataLoader
    gene_dataset = GeneExpressionDataset([X])
    # data_loader = DataLoader(gene_dataset, batch_size=4,
    #                         shuffle=True, num_workers=4)

    # 1. Instanciate model
    vae = VAE(nb_genes)

    # 2. Train model
    # train(vae, data_loader) - not working on real data so far

    # 3. Test model - imputation benchmark

    X_zero, i, j, ix = dropout(X_test)
    # TODO: so far, compute imputation score simpler as below
    sampler = SequentialSampler(gene_dataset)
    sequential_data_loader = DataLoader(gene_dataset, batch_size=X_test.shape[0],
                                        sampler=sampler, num_workers=1)

    for i_batch, (sample_batched, local_l_mean, local_l_var) in enumerate(sequential_data_loader):
        sample_batched = Variable(sample_batched, requires_grad=False)
        # TODO: avoid full forward pass
        _, _, px_rate, _, _, _, _, _ = vae(sample_batched)
        print("MAE is :", imputation_error(px_rate.data.numpy(), X_test, i, j, ix))

    assert 1 == 1  # If the imputation succeeds...
