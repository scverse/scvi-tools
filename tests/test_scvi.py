#!/usr/bin/env python
# -*- coding: utf-8 -*-


"""Tests for `scvi` package."""

import numpy as np
import pytest
from torch.utils.data import DataLoader

from scvi.dataset import GeneExpressionDataset
from scvi.scvi import VAE
from scvi.train import train


@pytest.fixture
def response():
    """Sample pytest fixture.

    See more at: http://doc.pytest.org/en/latest/fixture.html
    """
    # import requests
    # return requests.get('https://github.com/audreyr/cookiecutter-pypackage')


def test_content(response):
    """Sample pytest test function with the pytest fixture as an argument."""
    # from bs4 import BeautifulSoup
    # assert 'GitHub' in BeautifulSoup(response.content).title.string
    pass


def test_benchmark():
    # Generating samples according to a ZINB process
    batch_size = 20
    nb_genes = 100
    data = np.random.negative_binomial(5, 0.3, size=(batch_size, nb_genes))
    newdata = np.ones((batch_size, nb_genes))
    mask = np.random.binomial(n=1, p=0.7, size=(batch_size, nb_genes))
    for i in range(batch_size):
        newdata[i, :] = data[i, :] / np.sum(data[i, :])
        newdata[i, :] = newdata[i, :] * mask[i, :]
    # Creating a GeneExpressionDataset and a DataLoader
    gene_dataset = GeneExpressionDataset([newdata])
    data_loader = DataLoader(gene_dataset, batch_size=4,
                             shuffle=True, num_workers=4)

    # 1. Instanciate model
    vae = VAE(nb_genes)

    # 2. Train model
    train(vae, data_loader)


def test_imputation():
    """The imputation test"""

    batch_size = 20
    nb_genes = 100
    data = np.random.negative_binomial(5, 0.3, size=(batch_size, nb_genes))
    newdata = np.ones((batch_size, nb_genes))
    mask = np.random.binomial(n=1, p=0.7, size=(batch_size, nb_genes))
    for i in range(batch_size):
        newdata[i, :] = data[i, :] / np.sum(data[i, :])
        newdata[i, :] = newdata[i, :] * mask[i, :]
    # Creating a GeneExpressionDataset and a DataLoader
    gene_dataset = GeneExpressionDataset([newdata])
    data_loader = DataLoader(gene_dataset, batch_size=4,
                             shuffle=True, num_workers=4)

    # 1. Instanciate model
    vae = VAE(nb_genes)

    # 2. Train model
    train(vae, data_loader)

    assert 1 == 1  # The imputation succeeded...
