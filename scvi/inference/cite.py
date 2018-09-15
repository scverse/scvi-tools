import copy

import matplotlib.pyplot as plt
import torch

from scvi.inference import Posterior
from scvi.inference import Trainer
from . import UnsupervisedTrainer

plt.switch_backend('agg')

class CitePosterior(Posterior):

    def ll_umi(self, verbose=False):
        self.mode = 'umi'
        ll = self.compute_log_likelihood(self.model, self)
        if verbose:
            print("LL UMI: %.4f" % ll)
        return ll

    def ll_adt(self, verbose=False):
        self.mode = 'adt'
        ll = self.compute_log_likelihood(self.model, self)
        if verbose:
            print("LL ADT: %.4f" % ll)
        return ll

    def ll(self, verbose=False):
        self.mode = 'total'
        ll = self.compute_log_likelihood(self.model, self)
        if verbose:
            print("LL UMI: %.4f" % ll)
        return ll

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
            self.test_set.to_monitor = ['ll_umi', 'll_adt']

    def loss(self, tensors):
        sample_batch, local_l_mean, local_l_var, batch_index, _ = tensors
        reconst_loss_umi, reconst_loss_adt, kl_divergence = self.model(sample_batch, local_l_mean, local_l_var, batch_index)
        loss = torch.mean(reconst_loss_umi + reconst_loss_adt + self.kl_weight * kl_divergence)
        return loss
