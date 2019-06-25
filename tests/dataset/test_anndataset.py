from unittest import TestCase

import anndata
import numpy as np
import torch

from scvi.dataset import AnnDatasetFromAnnData, DownloadableAnnDataset, GeneExpressionDataset
from scvi.inference import UnsupervisedTrainer
from scvi.models import VAE

use_cuda = torch.cuda.is_available()


class TestAnnDatasetFromAnnData(TestCase):
    def test_init(self):
        data = np.random.randint(1, 5, size=(3, 7))
        ad = anndata.AnnData(data)
        dataset = AnnDatasetFromAnnData(ad)
        self.assertEqual(3, dataset.nb_cells)
        self.assertEqual(7, dataset.nb_genes)

    def test_train_one(self):
        data = np.random.randint(1, 5, size=(4, 7))
        ad = anndata.AnnData(data)
        dataset = AnnDatasetFromAnnData(ad)
        unsupervised_training_one_epoch(dataset)


class TestDownloadableAnnDataset(TestCase):
    def test_populate(self):
        data = np.random.randint(1, 5, size=(3, 7))
        dataset = DownloadableAnnDataset(delayed_populating=True)
        dataset.populate_from_data(data)
        self.assertEqual(3, dataset.nb_cells)
        self.assertEqual(7, dataset.nb_genes)

    def test_train_one(self):
        dataset = DownloadableAnnDataset("TM_droplet_mat.h5ad", save_path="tests/data")
        unsupervised_training_one_epoch(dataset)


def unsupervised_training_one_epoch(dataset: GeneExpressionDataset):
    vae = VAE(dataset.nb_genes, dataset.n_batches, dataset.n_labels)
    trainer = UnsupervisedTrainer(vae, dataset, train_size=0.5, use_cuda=use_cuda)
    trainer.train(n_epochs=1)
