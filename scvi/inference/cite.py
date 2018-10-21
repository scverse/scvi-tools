import copy

import matplotlib.pyplot as plt
import torch
from torch.distributions import Normal, Poisson, LogNormal
import numpy as np

from scvi.inference import Posterior
from scvi.inference import Trainer
from . import UnsupervisedTrainer

plt.switch_backend('agg')

class CitePosterior(Posterior):

    def ll_umi(self, verbose=False):
        self.mode = 'umi'
        ll = self.compute_log_likelihood(self.model)
        if verbose:
            print("LL UMI: %.4f" % ll)
        return ll

    def ll_adt(self, verbose=False):
        self.mode = 'adt'
        ll = self.compute_log_likelihood(self.model)
        if verbose:
            print("LL ADT: %.4f" % ll)
        return ll

    def ll(self, verbose=False):
        self.mode = 'total'
        ll = self.compute_log_likelihood(self.model)
        if verbose:
            print("LL UMI: %.4f" % ll)
        return ll

    def mse_adt(self, verbose=False):
        mse, _ = self.compute_adt_error(self.model)
        if verbose:
            print("MSE ADT: %.4f" % mse)
        return mse

    def mae_adt(self, verbose=False):
        _, mae = self.compute_adt_error(self.model)
        if verbose:
            print("MAE ADT: %.4f" % mae)
        return mae

    def compute_log_likelihood(self, vae, **kwargs):
        r""" Computes log p(x/z), which is the reconstruction error .
            Differs from the marginal log likelihood, but still gives good
            insights on the modeling of the data, and is fast to compute
        """
        # Iterate once over the posterior and computes the total log_likelihood
        log_lkl = 0
        for i_batch, tensors in enumerate(self):
            sample_batch, local_l_mean, local_l_var, batch_index, labels = tensors[:5]  # general fish case
            reconst_loss_umi, reconst_loss_adt, kl_divergence = vae(sample_batch, local_l_mean, local_l_var, batch_index=batch_index,
                                              y=labels, **kwargs)
            if self.mode == 'umi':
                log_lkl += torch.sum(reconst_loss_umi).item()
            elif self.mode == 'adt':
                log_lkl += torch.sum(reconst_loss_adt).item()
            else:
                reconst_loss = reconst_loss_umi + reconst_loss_adt
                log_lkl += torch.sum(reconst_loss).item()

        n_samples = len(self.indices)
        return log_lkl / n_samples

    def compute_adt_error(self, vae, mc_samples=100, **kwargs):
        r""" TODO
        """
        # Iterate once over the posterior and computes the total log_likelihood
        # mean squared error
        mse = 0
        # mean absolute error
        mae = 0
        n_samples = 0
        for i_batch, tensors in enumerate(self):
            sample_batch, local_l_mean, local_l_var, batch_index, labels = tensors[:5]  # general fish case
            n_samples += sample_batch.shape[0]
            for i in range(mc_samples):
                means = vae.get_sample_rate(sample_batch, mode='adt')
                var_disp = vae.get_sample_dispersion(sample_batch, mode='adt')

                # Parametrized as r, p
                if i == 0:
                    samples = NegativeBinomial(total_count=var_disp, probs=means / (means + var_disp)).sample()
                else:
                    samples += NegativeBinomial(total_count=var_disp, probs=means / (means + var_disp)).sample()

            # Take average sample
            samples /= mc_samples

            adt_ = sample_batch[:, vae.n_input_genes:]
            error = torch.sum((samples - adt_)**2).item()
            mse += error

            mae += torch.sum(torch.abs(samples - adt_)).item()

        n_samples = len(self.indices)
        mse /= (n_samples * vae.n_input_proteins)
        mae /= (n_samples * vae.n_input_proteins)
        return mse, mae

    def imputation_list(self, n_samples=1):
        original_list = []
        imputed_list = []
        batch_size = 10000  # self.data_loader_kwargs['batch_size'] // n_samples
        for tensors, corrupted_tensors in zip(self.uncorrupted().sequential(batch_size=batch_size),
                                              self.corrupted().sequential(batch_size=batch_size)):
            batch = tensors[0]
            actual_batch_size = batch.size(0)
            dropout_batch, _, _, batch_index, labels = corrupted_tensors
            px_rate_umi = self.model.get_sample_rate_umi(dropout_batch, batch_index=batch_index, y=labels, n_samples=n_samples)
            px_rate_adt = self.model.get_sample_rate_adt(dropout_batch, batch_index=batch_index, y=labels, n_samples=n_samples)
            # concatenate the imputed counts
            px_rate = torch.cat((px_rate_umi, px_rate_adt), 1)

            indices_dropout = torch.nonzero(batch - dropout_batch)
            i = indices_dropout[:, 0]
            j = indices_dropout[:, 1]

            batch = batch.unsqueeze(0).expand((n_samples, batch.size(0), batch.size(1)))
            original = np.array(batch[:, i, j].view(-1).cpu())
            imputed = np.array(px_rate[..., i, j].view(-1).cpu())

            cells_index = np.tile(np.array(i.cpu()), n_samples)

            original_list += [original[cells_index == i] for i in range(actual_batch_size)]
            imputed_list += [imputed[cells_index == i] for i in range(actual_batch_size)]
        return original_list, imputed_list

    def differential_expression_stats(self, M_sampling=100, mode='adt'):
        """
        Output average over statistics in a symmetric way (a against b)
        forget the sets if permutation is True
        :param vae: The generative vae and encoder network
        :param data_loader: a data loader for a particular dataset
        :param M_sampling: number of samples
        :return: A 1-d vector of statistics of size n_genes
        """
        px_scales = []
        all_labels = []
        batch_size = max(self.data_loader_kwargs['batch_size'] // M_sampling, 2)  # Reduce batch_size on GPU
        for tensors in self.update({"batch_size": batch_size}):
            sample_batch, _, _, batch_index, labels = tensors
            px_scales += [
                np.array((self.model.get_sample_scale(
                    sample_batch, batch_index=batch_index, y=labels, n_samples=M_sampling, mode=mode)
                         ).cpu())]

            # Align the sampling
            if M_sampling > 1:
                px_scales[-1] = (px_scales[-1].transpose((1, 0, 2))).reshape(-1, px_scales[-1].shape[-1])
            all_labels += [np.array((labels.repeat(1, M_sampling).view(-1, 1)).cpu())]

        px_scales = np.concatenate(px_scales)
        all_labels = np.concatenate(all_labels).ravel()  # this will be used as boolean

        return px_scales, all_labels


class CiteTrainer(UnsupervisedTrainer):
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
    default_metrics_to_monitor = ['ll_umi', 'll_adt']

    def __init__(self, model, gene_dataset, train_size=0.8, test_size=None, kl=None, **kwargs):
        super(CiteTrainer, self).__init__(model, gene_dataset, train_size=0.8, test_size=None, kl=None, **kwargs)
        self.kl = kl
        if type(self) is CiteTrainer:
            self.train_set, self.test_set = self.train_test(model, gene_dataset, train_size, test_size, type_class=CitePosterior)
            self.train_set.to_monitor = ['ll_umi', 'll_adt']
            self.test_set.to_monitor = ['ll_umi', 'll_adt', 'mse_adt', 'mae_adt']

    def loss(self, tensors):
        sample_batch, local_l_mean, local_l_var, batch_index, _ = tensors
        reconst_loss_umi, reconst_loss_adt, kl_divergence = self.model(sample_batch, local_l_mean, local_l_var, batch_index)
        loss = torch.mean(reconst_loss_umi + reconst_loss_adt + self.kl_weight * kl_divergence)
        return loss


import torch.nn.functional as F
from torch.distributions import constraints
from torch.distributions.distribution import Distribution
from torch.distributions.utils import broadcast_all, probs_to_logits, lazy_property, logits_to_probs

# CODE FROM https://pytorch.org/docs/master/_modules/torch/distributions/negative_binomial.html
class NegativeBinomial(Distribution):
    r"""
    Creates a Negative Binomial distribution, i.e. distribution
    of the number of independent identical Bernoulli trials
    needed before :attr:`total_count` failures are achieved. The probability
    of success of each Bernoulli trial is :attr:`probs`.

    Args:
        total_count (float or Tensor): non-negative number of negative Bernoulli
            trials to stop, although the distribution is still valid for real
            valued count
        probs (Tensor): Event probabilities of success in the half open interval [0, 1)
        logits (Tensor): Event log-odds for probabilities of success
    """

    def __init__(self, total_count, probs=None, logits=None, validate_args=None):
        if (probs is None) == (logits is None):
            raise ValueError("Either `probs` or `logits` must be specified, but not both.")
        if probs is not None:
            self.total_count, self.probs, = broadcast_all(total_count, probs)
            self.total_count = self.total_count.type_as(self.probs)
        else:
            self.total_count, self.logits, = broadcast_all(total_count, logits)
            self.total_count = self.total_count.type_as(self.logits)

        self._param = self.probs if probs is not None else self.logits
        batch_shape = self._param.size()
        super(NegativeBinomial, self).__init__(batch_shape, validate_args=validate_args)

    def expand(self, batch_shape, _instance=None):
        new = self._get_checked_instance(NegativeBinomial, _instance)
        batch_shape = torch.Size(batch_shape)
        new.total_count = self.total_count.expand(batch_shape)
        if 'probs' in self.__dict__:
            new.probs = self.probs.expand(batch_shape)
            new._param = new.probs
        else:
            new.logits = self.logits.expand(batch_shape)
            new._param = new.logits
        super(NegativeBinomial, new).__init__(batch_shape, validate_args=False)
        new._validate_args = self._validate_args
        return new


    def _new(self, *args, **kwargs):
        return self._param.new(*args, **kwargs)

    @property
    def mean(self):
        return self.total_count * torch.exp(self.logits)

    @property
    def variance(self):
        return self.mean / torch.sigmoid(-self.logits)

    @lazy_property
    def logits(self):
        return probs_to_logits(self.probs, is_binary=True)


    @lazy_property
    def probs(self):
        return logits_to_probs(self.logits, is_binary=True)


    @property
    def param_shape(self):
        return self._param.size()

    @lazy_property
    def _gamma(self):
        return torch.distributions.Gamma(concentration=self.total_count,
                                         rate=torch.exp(-self.logits))

    def sample(self, sample_shape=torch.Size()):
        with torch.no_grad():
            rate = self._gamma.sample(sample_shape=sample_shape)
            return torch.poisson(rate)


    def log_prob(self, value):
        if self._validate_args:
            self._validate_sample(value)

        log_unnormalized_prob = (self.total_count * F.logsigmoid(-self.logits) +
                                 value * F.logsigmoid(self.logits))

        log_normalization = (-torch.lgamma(self.total_count + value) + torch.lgamma(1. + value) +
                             torch.lgamma(self.total_count))

        return log_unnormalized_prob - log_normalization
