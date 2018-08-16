import copy
from itertools import cycle

import matplotlib.pyplot as plt
import numpy as np
import torch
from torch.nn import functional as F

from scvi.dataset import GeneExpressionDataset
from scvi.inference import Trainer
from scvi.inference.trainer import ToCudaDataLoader
from .posterior import PosteriorFish, plot_imputation, Posterior

plt.switch_backend('agg')


class UnsupervisedTrainer(Trainer):
    r"""The VariationalInference class for the unsupervised training of an autoencoder.

    Args:
        :model: A model instance from class ``VAE``, ``VAEC``, ``SVAEC``
        :gene_dataset: A gene_dataset instance like ``CortexDataset()``
        :train_size: The train size, either a float between 0 and 1 or and integer for the number of training samples
         to use Default: ``0.8``.
        :\**kwargs: Other keywords arguments from the general Inference class.

    Examples:
        >>> gene_dataset = CortexDataset()
        >>> vae = VAE(gene_dataset.nb_genes, n_batch=gene_dataset.n_batches * False,
        ... n_labels=gene_dataset.n_labels)

        >>> infer = VariationalInference(gene_dataset, vae, train_size=0.5)
        >>> infer.train(n_epochs=20, lr=1e-3)
    """
    default_metrics_to_monitor = ['ll']

    def __init__(self, model, gene_dataset, train_size=0.8, **kwargs):
        super(UnsupervisedTrainer, self).__init__(model, gene_dataset, **kwargs)
        self.kl = None
        self.sequential = self.create_posterior(name='sequential')
        self.train_set, self.test_set = self.train_test(model, gene_dataset, train_size)
        self.train_set.to_monitor = ['ll']
        self.posteriors_dict = {'train_set': self.train_set, 'test_set': self.test_set}
        self.posteriors_loop = ['train_set']
        self.train_set.to_monitor = ['ll']
        self.test_set.to_monitor = ['ll']

    def loss(self, tensors):
        sample_batch, local_l_mean, local_l_var, batch_index, _ = tensors
        reconst_loss, kl_divergence = self.model(sample_batch, local_l_mean, local_l_var, batch_index)
        loss = torch.mean(reconst_loss + self.kl_weight * kl_divergence)
        return loss

    def on_epoch_begin(self):
        self.kl_weight = self.kl if self.kl is not None else min(1, self.epoch / self.n_epochs)


class AdapterTrainer(UnsupervisedTrainer):
    def __init__(self, model, gene_dataset, posterior_test, frequency=5):
        super(AdapterTrainer, self).__init__(model, gene_dataset, frequency=frequency)
        self.test_set = posterior_test
        self.test_set.name = 'test_set'
        self.test_set.to_monitor = ['ll']
        self.posteriors_dict = {'test_set': self.test_set}
        self.params = list(self.model.z_encoder.parameters()) + list(self.model.l_encoder.parameters())
        self.z_encoder_state = copy.deepcopy(model.z_encoder.state_dict())
        self.l_encoder_state = copy.deepcopy(model.l_encoder.state_dict())

    def data_loaders_loop(self):
        return zip(self.test_set.data_loader)

    def train(self, n_path=10, n_epochs=50, **kwargs):
        # Training the model

        for i in range(n_path):
            # Re-initialize to create new path
            self.model.z_encoder.load_state_dict(self.z_encoder_state)
            self.model.l_encoder.load_state_dict(self.l_encoder_state)
            super(AdapterTrainer, self).train(n_epochs, params=self.params, **kwargs)

        return min(self.history["ll_test_set"])


class ImputationTrainer:
    def __init__(self, trainer, gene_dataset, rate=0.1, corruption="uniform"):
        self.corruption = corruption
        self.rate = rate
        corrupted_data = copy.deepcopy(gene_dataset.X)

        if corruption == "uniform":  # multiply the entry n with a Ber(0.9) random variable.
            i, j = np.nonzero(corrupted_data)
            ix = np.random.choice(range(len(i)), int(np.floor(rate * len(i))), replace=False)
            i, j = i[ix], j[ix]
            corrupted_data[i, j] *= np.random.binomial(n=np.ones(len(ix), dtype=np.int64), p=0.9)
        elif corruption == "binomial":  # multiply the entry n with a Bin(n, 0.9) random variable.
            i, j = (k.ravel() for k in np.indices(corrupted_data.shape))
            ix = np.random.choice(range(len(i)), int(np.floor(rate * len(i))), replace=False)
            i, j = i[ix], j[ix]
            corrupted_data[i, j] = np.random.binomial(n=corrupted_data[i, j].astype(np.int64), p=0.2)

        gene_dataset = GeneExpressionDataset(
            *GeneExpressionDataset.get_attributes_from_matrix(
                corrupted_data,
                batch_indices=gene_dataset.batch_indices,
                labels=gene_dataset.labels
            )
        )
        self.trainer = trainer
        self.posteriors_dict = trainer.posteriors_dict
        self.corrupted_posteriors_dict = dict()
        for key, posterior in self.posteriors_dict.items():
            kwargs = copy.copy(self.trainer.data_loaders_kwargs)
            kwargs['collate_fn'] = gene_dataset.collate_fn
            kwargs['sampler'] = copy.copy(posterior.data_loader.sampler)
            self.corrupted_posteriors_dict[key] = Posterior(trainer.model, gene_dataset,
                                                            ToCudaDataLoader(gene_dataset, **kwargs),
                                                            name=key)
        self.trainer.posteriors_dict = self.corrupted_posteriors_dict
        self.trainer.posteriors_loop = ['test_set']

    def train(self, *args, **kwargs):
        self.trainer.train(*args, **kwargs)

    def imputation_task(self, name, n_samples=1):
        original_list = []
        imputed_list = []
        batch_size = self.posteriors_dict[name].data_loader.batch_size // n_samples
        for tensors, corrupted_tensors in \
            zip((self.posteriors_dict[name].data_loader).sequential(batch_size=batch_size),
                (self.corrupted_posteriors_dict[name].data_loader).sequential(batch_size=batch_size)):
            batch = tensors[0]
            actual_batch_size = batch.size(0)
            dropout_batch, _, _, batch_index, labels = corrupted_tensors
            px_rate = self.trainer.model.get_sample_rate(dropout_batch, batch_index=batch_index, y=labels,
                                                         n_samples=n_samples)

            indices_dropout = torch.nonzero(batch - dropout_batch)
            i = indices_dropout[:, 0]
            j = indices_dropout[:, 1]

            batch = batch.unsqueeze(0).expand((n_samples, batch.size(0), batch.size(1)))
            original = np.array(batch[:, i, j].view(-1).cpu())
            imputed = np.array(px_rate[:, i, j].view(-1).cpu())

            cells_index = np.tile(np.array(i.cpu()), n_samples)

            original_list += [original[cells_index == i] for i in range(actual_batch_size)]
            imputed_list += [imputed[cells_index == i] for i in range(actual_batch_size)]

        return original_list, imputed_list

    def imputation(self, name, verbose=False, n_samples=1, title="Imputation"):
        original_list, imputed_list = self.imputation_task(name, n_samples=n_samples)
        # Median of medians for all distances
        median_imputation_score = np.median(np.abs(np.concatenate(original_list) - np.concatenate(imputed_list)))

        # Mean of medians for each cell
        imputation_cells = []
        for original, imputed in zip(original_list, imputed_list):
            has_imputation = len(original) and len(imputed)
            imputation_cells += [np.median(np.abs(original - imputed)) if has_imputation else 0]
        mean_imputation_score = np.mean(imputation_cells)

        if verbose:
            print("Imputation Scores [corruption:%s - rate:%.2f] on  %s after %i:"
                  "\nMedian of Median: %.4f\nMean of Median for each cell: %.4f" %
                  (self.corruption, self.rate, name, self.n_epochs, median_imputation_score, mean_imputation_score))

        plot_imputation(np.concatenate(original_list), np.concatenate(imputed_list), title=title)
        return original_list, imputed_list


class TrainerFish(Trainer):
    r"""The VariationalInference class for the unsupervised training of an autoencoder.

    Args:
        :model: A model instance from class ``VAE``, ``VAEC``, ``SVAEC``
        :gene_dataset: A gene_dataset instance like ``CortexDataset()``
        :train_size: The train size, either a float between 0 and 1 or and integer for the number of training samples
         to use Default: ``0.8``.
        :**kwargs: Other keywords arguments from the general Inference class.

    Examples:
        >>> gene_dataset_seq = CortexDataset()
        >>> gene_dataset_fish = SmfishDataset()
        >>> vae = VAE(gene_dataset_seq.nb_genes, gene_dataset_fish.nb_genes,
        ... n_labels=gene_dataset.n_labels, use_cuda=True)

        >>> infer = VariationalInference(gene_dataset_seq, gene_dataset_fish, vae, train_size=0.5)
        >>> infer.train(n_epochs=20, lr=1e-3)
    """
    default_metrics_to_monitor = ['ll']

    def __init__(self, model, gene_dataset_seq, gene_dataset_fish, train_size=0.8, test_size=None,
                 use_cuda=True, cl_ratio=0, n_epochs_even=1, n_epochs_kl=2000, n_epochs_cl=1, seed=0, **kwargs):
        super(TrainerFish, self).__init__(model, gene_dataset_seq, use_cuda=use_cuda, **kwargs)
        self.kl = None
        self.cl_ratio = cl_ratio
        self.n_epochs_cl = n_epochs_cl
        self.n_epochs_even = n_epochs_even
        self.n_epochs_kl = n_epochs_kl
        self.weighting = 0
        self.kl_weight = 0
        self.classification_ponderation = 0

        self.train_seq, self.test_seq = self.train_test(self.model, gene_dataset_seq, train_size, test_size, seed,
                                                        suffix='_seq')
        self.train_fish, self.test_fish = self.train_test(self.model, gene_dataset_fish, train_size, test_size, seed,
                                                          suffix='_fish', cls=PosteriorFish)
        self.test_seq.to_monitor = ['ll']
        self.test_fish.to_monitor = ['ll']
        self.posteriors = [self.test_seq, self.test_fish]

    def data_loaders_loop(self):
        posteriors_loop = [self.train_seq, self.train_fish]
        data_loaders_loop = [posterior.data_loader for posterior in posteriors_loop]
        return zip(data_loaders_loop[0], *[cycle(data_loader) for data_loader in data_loaders_loop[1:]])

    def loss(self, tensors_seq, tensors_fish):
        sample_batch, local_l_mean, local_l_var, batch_index, labels = tensors_seq
        reconst_loss, kl_divergence = self.model(sample_batch, local_l_mean, local_l_var, batch_index, mode="scRNA",
                                                 weighting=self.weighting)
        # If we want to add a classification loss
        if self.cl_ratio != 0:
            reconst_loss += self.cl_ratio * F.cross_entropy(self.model.classify(sample_batch, mode="scRNA"),
                                                            labels.view(-1))
        loss = torch.mean(reconst_loss + self.kl_weight * kl_divergence)
        sample_batch_fish, local_l_mean, local_l_var, batch_index_fish, _, _, _ = tensors_fish
        reconst_loss_fish, kl_divergence_fish = self.model(sample_batch_fish, local_l_mean, local_l_var,
                                                           batch_index_fish, mode="smFISH")
        loss_fish = torch.mean(reconst_loss_fish + self.kl_weight * kl_divergence_fish)
        loss = loss * sample_batch.size(0) + loss_fish * sample_batch_fish.size(0)
        loss /= (sample_batch.size(0) + sample_batch_fish.size(0))
        return loss + loss_fish

    def on_epoch_begin(self):
        self.weighting = min(1, self.epoch / self.n_epochs_even)
        self.kl_weight = self.kl if self.kl is not None else min(1, self.epoch / self.n_epochs_kl)
        self.classification_ponderation = min(1, self.epoch / self.n_epochs_cl)
