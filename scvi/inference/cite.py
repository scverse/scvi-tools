import matplotlib.pyplot as plt
import torch
from torch.distributions import Poisson, LogNormal, Gamma
import numpy as np

from scvi.inference import Posterior
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
            print("LL UMI+ADT: %.4f" % ll)
        return ll

    def mse_adt(self, verbose=False):
        mse, _, _ = self.compute_adt_error(self.model)
        if verbose:
            print("MSE ADT: %.4f" % mse)
        return mse

    def mae_adt(self, verbose=False):
        _, mae, _ = self.compute_adt_error(self.model)
        if verbose:
            print("MAE ADT: %.4f" % mae)
        return mae

    def compute_log_likelihood(self, vae, **kwargs):
        r""" Computes log p(x/z), which is the reconstruction error .
            Differs from the marginal log likelihood, but still gives good
            insights on the modeling of the data, and is fast to compute

            This is really a helper function to self.ll, self.ll_adt, etc.
        """
        # Iterate once over the posterior and computes the total log_likelihood
        log_lkl = 0
        for i_batch, tensors in enumerate(self):
            # general fish case
            sample_batch, local_l_mean, local_l_var, batch_index, labels = tensors[:5]
            reconst_loss_umi, reconst_loss_adt, kl_divergence = vae(sample_batch, local_l_mean,
                                                                    local_l_var, batch_index=batch_index,
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

    def compute_adt_error(self, vae, mc_samples=1, **kwargs):
        r""" Samples from the posterior predictive p(x_new|x) mc_samples times
             and takes the average. So, when mc_samples=1, it's one sample from
             the posterior predictive. The given 'x' represents the cells in this
             posterior class.

             Returns the MSE, MAE between new sample test data, as well as LL of
             new sample test data.

        """
        # Iterate once over the posterior and computes the total log_likelihood
        # mean squared error
        mse = 0
        # mean absolute error
        mae = 0
        # log likelihood
        log_lkl = 0
        for i_batch, tensors in enumerate(self):
            # general fish case
            sample_batch, local_l_mean, local_l_var, batch_index, labels = tensors[:5]
            # For loop could probably be removed through n_samples parameter
            for i in range(mc_samples):
                means = vae.get_sample_rate(sample_batch, mode='adt')
                var_disp = vae.get_sample_dispersion(sample_batch, mode='adt')

                # Parametrized as r, p
                if vae.reconstruction_loss_adt == 'nb':
                    r = var_disp
                    p = means / (means + var_disp)
                    odds = p / (1 - p)
                    rate = Gamma(concentration=r, rate=1 / odds).sample()
                    if i == 0:
                        # samples = NegativeBinomial(total_count=var_disp, probs=means / (means + var_disp)).sample()
                        samples = torch.poisson(rate)
                    else:
                        # samples += NegativeBinomial(total_count=var_disp, probs=means / (means + var_disp)).sample()
                        samples += torch.poisson(rate)
                elif vae.reconstruction_loss_adt == 'poisson':
                    if i == 0:
                        samples = Poisson(means).sample()
                    else:
                        samples += Poisson(means).sample()
                # log_normal
                else:
                    if i == 0:
                        samples = LogNormal(means, var_disp).sample()
                    else:
                        samples += LogNormal(means, var_disp).sample()

            # Take average sample
            samples /= mc_samples

            umi_ = sample_batch[:, :vae.n_input_genes]
            reconst_loss_umi, reconst_loss_adt, kl_divergence = vae(torch.cat((umi_, samples), 1), local_l_mean,
                                                                    local_l_var, batch_index=batch_index,
                                                                    y=labels, **kwargs)
            log_lkl += torch.sum(reconst_loss_adt).item()

            adt_ = sample_batch[:, vae.n_input_genes:]
            error = torch.sum((samples - adt_)**2, dim=0)
            mse += error

            mae += torch.sum(torch.abs(samples - adt_), dim=0)

        n_samples = len(self.indices)
        mse /= n_samples
        mae /= n_samples
        log_lkl /= n_samples
        return mse, mae, log_lkl

    def differential_expression_stats(self, M_sampling=100, mode='adt'):
        r"""
        Output average over statistics in a symmetric way (a against b)
        forget the sets if permutation is True
        :param vae: The generative vae and encoder network
        :param data_loader: a data loader for a particular dataset
        :param M_sampling: number of samples
        :param mode: UMI or ADT data
        :return: A 1-d vector of statistics of size n_genes
        """
        px_scales = []
        all_labels = []
        # Reduce batch_size on GPU
        batch_size = max(
            self.data_loader_kwargs['batch_size'] // M_sampling, 2)
        for tensors in self.update({"batch_size": batch_size}):
            sample_batch, _, _, batch_index, labels = tensors
            px_scales += [
                np.array((self.model.get_sample_scale(
                    sample_batch, batch_index=batch_index, y=labels, n_samples=M_sampling, mode=mode)
                ).cpu())]

            # Align the sampling
            if M_sampling > 1:
                px_scales[-1] = (px_scales[-1].transpose((1, 0, 2))
                                 ).reshape(-1, px_scales[-1].shape[-1])
            all_labels += [np.array((labels.repeat(1,
                                                   M_sampling).view(-1, 1)).cpu())]

        px_scales = np.concatenate(px_scales)
        # this will be used as boolean
        all_labels = np.concatenate(all_labels).ravel()

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
        super(CiteTrainer, self).__init__(model, gene_dataset,
                                          train_size=0.8, test_size=None, kl=None, **kwargs)
        self.kl = kl
        if type(self) is CiteTrainer:
            self.train_set, self.test_set = self.train_test(
                model, gene_dataset, train_size, test_size, type_class=CitePosterior)
            self.train_set.to_monitor = ['ll_umi', 'll_adt']
            self.test_set.to_monitor = ['ll_umi', 'll_adt']

    def loss(self, tensors):
        sample_batch, local_l_mean, local_l_var, batch_index, _ = tensors
        reconst_loss_umi, reconst_loss_adt, kl_divergence = self.model(
            sample_batch, local_l_mean, local_l_var, batch_index)
        loss = torch.mean(reconst_loss_umi + reconst_loss_adt +
                          self.kl_weight * kl_divergence)
        return loss
