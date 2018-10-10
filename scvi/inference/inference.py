import copy

import matplotlib.pyplot as plt
import torch
from scvi.models.classifier import Classifier
import torch.nn.functional as F


from . import Trainer

plt.switch_backend('agg')


class UnsupervisedTrainer(Trainer):
    r"""The VariationalInference class for the unsupervised training of an autoencoder.

    Args:
        :model: A model instance from class ``VAE``, ``VAEC``, ``SCANVI``
        :gene_dataset: A gene_dataset instance like ``CortexDataset()``
        :train_size: The train size, either a float between 0 and 1 or and integer for the number of training samples
         to use Default: ``0.8``.
        :\*\*kwargs: Other keywords arguments from the general Trainer class.

    Examples:
        >>> gene_dataset = CortexDataset()
        >>> vae = VAE(gene_dataset.nb_genes, n_batch=gene_dataset.n_batches * False,
        ... n_labels=gene_dataset.n_labels)

        >>> infer = VariationalInference(gene_dataset, vae, train_size=0.5)
        >>> infer.train(n_epochs=20, lr=1e-3)
    """
    default_metrics_to_monitor = ['ll']

    def __init__(self, model, gene_dataset, train_size=0.8, test_size=None, kl=None, **kwargs):
        super(UnsupervisedTrainer, self).__init__(model, gene_dataset, **kwargs)
        self.kl = kl
        if type(self) is UnsupervisedTrainer:
            self.train_set, self.test_set = self.train_test(model, gene_dataset, train_size, test_size)
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
        self.kl_weight = self.kl if self.kl is not None else min(1, self.epoch / 400)  # self.n_epochs)


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


class AdversarialTrainerVAE(Trainer):
    r"""The modified UnsupervisedTrainer class for the unsupervised training of an autoencoder.
    """
    default_metrics_to_monitor = ['ll']

    def __init__(self, model, gene_dataset, train_size=0.8, test_size=None,
                 n_epochs_even=1, n_epochs_cl=1, warm_up=10,
                 scale=50, **kwargs):
        super(AdversarialTrainerVAE, self).__init__(model, gene_dataset, **kwargs)
        print("I am the adversarial Trainer")
        self.kl = None
        self.n_epochs_cl = n_epochs_cl
        self.n_epochs_even = n_epochs_even
        self.weighting = 0
        self.kl_weight = 0
        self.classification_ponderation = 0
        self.warm_up = warm_up
        self.scale = scale

        self.train_set, self.test_set = self.train_test(model, gene_dataset, train_size, test_size)
        self.train_set.to_monitor = ['ll']
        self.test_set.to_monitor = ['ll']


    @property
    def posteriors_loop(self):
        return ['train_set']

    def train(self, n_epochs=20, lr=1e-3, weight_decay=1e-6, params=None):
        self.adversarial_cls = Classifier(self.model.n_latent, n_labels=self.model.n_batch, n_layers=3)
        if self.use_cuda:
            self.adversarial_cls.cuda()
        self.optimizer_cls = torch.optim.Adam(filter(lambda p: p.requires_grad, self.adversarial_cls.parameters()),
                                              lr=lr, weight_decay=weight_decay)
        super(AdversarialTrainerVAE, self).train(n_epochs=n_epochs, lr=lr, params=None)

    def loss(self, tensors):
        sample_batch, local_l_mean, local_l_var, batch_index, _ = tensors
        reconst_loss, kl_divergence = self.model(sample_batch, local_l_mean, local_l_var, batch_index)
        loss = torch.mean(reconst_loss + self.kl_weight * kl_divergence)
        if self.epoch > self.warm_up:
            z = self.model.sample_from_posterior_z(sample_batch)
            cls_loss = (self.scale * F.cross_entropy(self.adversarial_cls(z), torch.zeros_like(batch_index).view(-1)))
            self.optimizer_cls.zero_grad()
            cls_loss.backward(retain_graph=True)
            self.optimizer_cls.step()
        else:
            cls_loss = 0
        return loss - cls_loss

    def on_epoch_begin(self):
        self.kl_weight = self.kl if self.kl is not None else min(1, self.epoch / 400)

