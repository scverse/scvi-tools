import copy

import matplotlib.pyplot as plt
import torch
from torch.nn import functional as F

from scvi.inference import Trainer

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

    def __init__(self, model, gene_dataset, train_size=0.8, kl=None, **kwargs):
        super(UnsupervisedTrainer, self).__init__(model, gene_dataset, **kwargs)
        self.kl = kl
        if type(self) is UnsupervisedTrainer:
            self.train_set, self.test_set = self.train_test(model, gene_dataset, train_size)
            self.train_set.to_monitor = ['ll']
            self.test_set.to_monitor = ['ll']

    @property
    def posteriors_loop(self):
        return ['train_set']

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
        self.test_set.to_monitor = ['ll']
        self.params = list(self.model.z_encoder.parameters()) + list(self.model.l_encoder.parameters())
        self.z_encoder_state = copy.deepcopy(model.z_encoder.state_dict())
        self.l_encoder_state = copy.deepcopy(model.l_encoder.state_dict())

    @property
    def posteriors_loop(self):
        return ['test_set']

    def train(self, n_path=10, n_epochs=50, **kwargs):
        for i in range(n_path):
            # Re-initialize to create new path
            self.model.z_encoder.load_state_dict(self.z_encoder_state)
            self.model.l_encoder.load_state_dict(self.l_encoder_state)
            super(AdapterTrainer, self).train(n_epochs, params=self.params, **kwargs)

        return min(self.history["ll_test_set"])


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

        self.train_seq, self.test_seq = self.train_test(self.model, gene_dataset_seq, train_size, test_size, seed)
        self.train_fish, self.test_fish = self.train_test(self.model, gene_dataset_fish, train_size, test_size, seed)
        self.test_seq.to_monitor = ['ll']
        self.test_fish.to_monitor = ['ll_fish']

    @property
    def posteriors_loop(self):
        return ['train_seq', 'train_fish']

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
