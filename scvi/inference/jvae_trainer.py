import logging
import sys
import time
from collections import defaultdict
from itertools import cycle
from functools import partial

import numpy as np
import torch
from torch import nn
from tqdm import trange

from scvi.inference import Posterior
from scvi.inference import Trainer
from scvi.models.log_likelihood import compute_elbo

logger = logging.getLogger(__name__)


class JPosterior(Posterior):
    def __init__(self, *args, mode=0, **kwargs):
        super().__init__(*args, **kwargs)
        self.mode = mode

    def elbo(self, verbose: bool=False) -> float:
        elbo = compute_elbo(self.model, self, mode=self.mode)
        if verbose:
            logger.info("ELBO : %.4f" % elbo)
        return elbo


class JVAETrainer(Trainer):
    r"""
    The VariationalInference class for the unsupervised training of an autoencoder.

        Args:
            :param model: A model instance from class ``JVAE``
            :param discriminator: A model instance of a classifier (with logit output)
            :param gene_dataset_list: list of gene_dataset instance like ``[CortexDataset(), SmfishDataset()]``
            :param train_size: Train size on cells
            :param kappa: float to weight the discriminator loss
            :param n_epochs_kl_warmup: Number of epochs for linear warmup of KL(q(z|x)||p(z)) term. After `n_epochs_kl_warmup`,
                the training objective is the ELBO. This might be used to prevent inactivity of latent units, and/or to
                improve clustering of latent space, as a long warmup turns the model into something more of an autoencoder.
            :param kwargs: Other keywords arguments from the general Trainer class.
    """

    default_metrics_to_monitor = ["elbo"]

    def __init__(
        self,
        model: nn.Module,
        discriminator: nn.Module,
        gene_dataset_list,
        train_size: float=0.8,
        use_cuda: bool=True,
        kappa: float=1.0,
        n_epochs_kl_warmup: int=400,
        **kwargs
    ):

        super().__init__(model, gene_dataset_list[0], use_cuda=use_cuda, **kwargs)
        self.n_epochs_kl_warmup = n_epochs_kl_warmup
        self.kappa_target = kappa
        self.all_dataset = [
            self.create_posterior(gene_dataset=gd, type_class=partial(JPosterior, mode=i))
            for i, gd in enumerate(gene_dataset_list)
        ]
        self.n_dataset = len(self.all_dataset)
        self.all_train, self.all_test = list(
            zip(*[self.train_test(model, gd, train_size, type_class=partial(JPosterior, mode=i))
                  for i, gd in enumerate(gene_dataset_list)])
        )
        for i, d in enumerate(self.all_train):
            self.register_posterior("train_%d" % i, d)
            d.to_monitor = ["elbo"]

        for i, d in enumerate(self.all_test):
            self.register_posterior("test_%d" % i, d)
            d.to_monitor = ['elbo']

        self.discriminator = discriminator
        if self.use_cuda:
            self.discriminator.cuda()

        self.kl_weight = None
        self.kappa = None
        self.compute_metrics_time = None
        self.n_epochs = None

        self.track_disc = []

    def loss_discriminator(self, tensors, predict_true_class=True):
        log_softmax = nn.LogSoftmax(dim=1)
        n_classes = self.n_dataset
        losses = []
        for i, z in enumerate(tensors):
            cls_logits = log_softmax(self.discriminator(z))

            if predict_true_class:
                cls_target = torch.zeros(n_classes, dtype=torch.float32, device=z.device)
                cls_target[i] = 1.0
            else:
                cls_target = torch.ones(n_classes, dtype=torch.float32, device=z.device) / (n_classes - 1)
                cls_target[i] = 0.0

            l_soft = cls_logits * cls_target
            cls_loss = -l_soft.sum(dim=1).mean()
            losses.append(cls_loss)

        total_loss = torch.stack(losses).sum()
        return total_loss

    def train(self, n_epochs=20, lr_d=1e-3, lr_g=1e-3, eps=0.01):
        self.compute_metrics_time = 0
        self.n_epochs = n_epochs

        self.compute_metrics()
        begin = time.time()

        with torch.set_grad_enabled(True):
            self.model.train()
            self.discriminator.train()

            d_params = filter(lambda p: p.requires_grad, self.discriminator.parameters())
            d_optimizer = torch.optim.Adam(d_params, lr=lr_d, eps=eps)

            g_params = filter(lambda p: p.requires_grad, self.model.parameters())
            g_optimizer = torch.optim.Adam(g_params, lr=lr_g, eps=eps)

            train_discriminator = self.n_dataset > 1 and self.kappa_target > 0

            with trange(n_epochs, desc="training", file=sys.stdout, disable=self.verbose) as progress_bar:
                for self.epoch in progress_bar:
                    self.on_epoch_begin()
                    progress_bar.update(1)
                    for tensors in self.data_loaders_loop():
                        if train_discriminator:
                            latent_tensors = []
                            for (i, (sample_batch, *_)) in enumerate(tensors):
                                z = self.model.sample_from_posterior_z(sample_batch, mode=i, deterministic=True)
                                latent_tensors.append(z)

                            # Train discriminator
                            d_loss = self.loss_discriminator([t.detach() for t in latent_tensors], True)
                            d_loss *= self.kappa
                            d_optimizer.zero_grad()
                            d_loss.backward()
                            d_optimizer.step()

                            # Train generative model to fool discriminator
                            fool_loss = self.loss_discriminator(latent_tensors, False)
                            fool_loss *= self.kappa
                            g_optimizer.zero_grad()
                            fool_loss.backward()
                            g_optimizer.step()

                        # Train generative model
                        g_loss = self.loss(tensors)
                        g_optimizer.zero_grad()
                        g_loss.backward()
                        g_optimizer.step()

                    if not self.on_epoch_end():
                        break

        self.model.eval()
        self.training_time += (time.time() - begin) - self.compute_metrics_time
        if self.verbose and self.frequency:
            print("\nTraining time:  %i s. / %i epochs" % (int(self.training_time), self.n_epochs))

    @property
    def posteriors_loop(self):
        return ['train_%d' % i for i in range(self.n_dataset)]

    def data_loaders_loop(self):
        posteriors = [self._posteriors[name] for name in self.posteriors_loop]
        # find the largest dataset to cycle over the others
        largest = np.argmax([posterior.gene_dataset.X.shape[0] for posterior in posteriors])

        data_loaders = [
            posterior if i == largest else cycle(posterior) for i, posterior in enumerate(posteriors)
        ]

        return zip(*data_loaders)

    def loss(self, tensors, return_details=False):
        reconstruction_losses = []
        kl_divergences = []
        losses = []
        total_batch_size = 0
        for (i, (sample_batch, local_l_mean, local_l_var, batch_index, labels, *_)) in enumerate(tensors):
            reconstruction_loss, kl_divergence = self.model(
                sample_batch, local_l_mean, local_l_var, batch_index, mode=i
            )
            loss = torch.mean(reconstruction_loss + self.kl_weight * kl_divergence) * sample_batch.size(0)
            total_batch_size += sample_batch.size(0)
            losses.append(loss)
            reconstruction_losses.append(reconstruction_loss)
            kl_divergences.append(kl_divergence)

        if return_details:
            return reconstruction_losses, kl_divergences

        averaged_loss = torch.stack(losses).sum() / total_batch_size
        return averaged_loss

    def compute_accuracy(self):
        confusion = []
        for i, posterior in enumerate(self.all_dataset):
            data = torch.from_numpy(posterior.gene_dataset.X)
            if self.use_cuda:
                data.cuda()

            z = self.model.sample_from_posterior_z(data, mode=i, deterministic=True)
            cls_z = nn.Softmax(dim=1)(self.discriminator(z)).detach()

            if self.use_cuda:
                cls_z = cls_z.cpu()
            cls_z = cls_z.numpy()

            row = cls_z.mean(axis=0)
            confusion.append(row)
        return np.array(confusion)

    def get_loss_magnitude(self, one_sample=True):
        total_reconstruction = defaultdict(float)
        total_kl_divergence = defaultdict(float)
        for tensors_list in self.data_loaders_loop():
            reconstruction_losses, kl_divergences = self.loss(tensors_list, return_details=True)
            for i in range(len(reconstruction_losses)):
                total_reconstruction[i] += float(reconstruction_losses[i].detach().numpy())
            for i in range(len(kl_divergences)):
                total_kl_divergence[i] += float(kl_divergences[i].detach().numpy())
            if one_sample:
                break
        return total_reconstruction, total_kl_divergence

    def on_epoch_begin(self):
        if self.n_epochs_kl_warmup is not None:
            self.kl_weight = min(1, self.epoch / self.n_epochs_kl_warmup)
        else:
            self.kl_weight = 1.0
        self.kappa = (self.epoch / self.n_epochs) * self.kappa_target

    def get_latent(self, deterministic=True):
        self.model.eval()
        latents = []
        for mode, dataset in enumerate(self.all_dataset):
            latent = []
            for tensors in dataset:
                sample_batch, local_l_mean, local_l_var, batch_index, label, *_ = tensors
                latent.append(self.model.sample_from_posterior_z(sample_batch, mode, deterministic=deterministic))

            latent = torch.cat(latent).detach().numpy()
            latents.append(latent)

        return latents

    def get_imputed_values(self, deterministic=True, normalized=True, decode_mode=None):
        self.model.eval()
        imputed_values = []
        for mode, dataset in enumerate(self.all_dataset):
            imputed_value = []
            for tensors in dataset:
                sample_batch, local_l_mean, local_l_var, batch_index, label, *_ = tensors
                if normalized:
                    imputed_value.append(
                        self.model.sample_scale(sample_batch, mode, batch_index, label, deterministic=deterministic,
                                                decode_mode=decode_mode)
                    )
                else:
                    imputed_value.append(
                        self.model.sample_rate(sample_batch, mode, batch_index, label, deterministic=deterministic,
                                               decode_mode=decode_mode)
                    )

            imputed_value = torch.cat(imputed_value).detach().numpy()
            imputed_values.append(imputed_value)

        return imputed_values
