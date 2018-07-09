from abc import abstractmethod

import numpy as np
import torch

from scvi.dataset import CortexDataset
from scvi.metrics.clustering import get_latent, entropy_batch_mixing
from .classification import compute_accuracy
from .differential_expression import de_stats, de_cortex
from .imputation import imputation
from .log_likelihood import compute_log_likelihood
from .visualization import show_t_sne

__all__ = ['compute_log_likelihood',
           'compute_accuracy',
           'show_t_sne']


class Metric:
    def __init__(self, mode):
        self.mode = mode

    @abstractmethod
    def __call__(self, infer, data_loader, name, *args, **kwargs):
        pass


class LLMetric(Metric):
    def __init__(self):
        super(LLMetric, self).__init__(mode='min')

    def __call__(self, infer, data_loader, name, *args, **kwargs):
        ll = compute_log_likelihood(infer.model, data_loader, *args, **kwargs)
        print("LL for %s is : %.4f" % (name, ll))
        return ll


class ImputationMetric(Metric):
    def __init__(self):
        self.imputation_task = ImputationTask()
        super(ImputationMetric, self).__init__(mode='min')

    def __call__(self, infer, data_loader, name, *args, **kwargs):
        imputation_score = torch.median(self.imputation_task(infer, data_loader, name, *args, **kwargs)).item()
        print("Median Imputation Score for %s is : %.4f" % (name, imputation_score))
        return imputation_score


class DEMetric(Metric):
    def __init__(self):
        self.de_task = DEStatsTask()
        super(DEMetric, self).__init__(mode='max')

    def __call__(self, infer, data_loader, name, *args, **kwargs):
        px_scale, all_labels = self.de_task(infer, data_loader, name, *args, **kwargs)

        if type(infer.gene_dataset) == CortexDataset:
            if 'use_cuda' in kwargs:
                kwargs.pop('use_cuda')
            de_score = de_cortex(px_scale, all_labels, infer.gene_dataset.gene_names, **kwargs)
            print("DE score for cortex on %s is : %.4f" % (name, de_score))
        return de_score


class BEMetric(Metric):
    def __init__(self):
        super(BEMetric, self).__init__(mode='max')

    def __call__(self, infer, data_loader, name, *args, **kwargs):
        if infer.gene_dataset.n_batches == 2:
            latent, batch_indices, labels = get_latent(infer.model, data_loader, *args, **kwargs)
            be_score = entropy_batch_mixing(latent, batch_indices)
            print("Entropy batch mixing :", be_score)
            return be_score


class AccuracyMetric(Metric):
    def __init__(self):
        super(AccuracyMetric, self).__init__(mode='max')

    def __call__(self, infer, data_loader, name, *args, **kwargs):
        acc = compute_accuracy(infer, data_loader, *args, **kwargs)
        print("Acc for %s is : %.4f" % (name, acc))
        return acc


class Task:
    @abstractmethod
    def __call__(self, infer, data_loader, name, *args, **kwargs):
        pass


class ImputationTask(Task):
    def __call__(self, infer, data_loader, name, *args, **kwargs):
        return imputation(infer.model, data_loader, *args, **kwargs)


class DEStatsTask(Task):
    def __call__(self, infer, data_loader, name, *args, **kwargs):
        return de_stats(infer.model, data_loader, *args, **kwargs)


class TsneTask(Task):
    def __call__(self, infer, data_loader, name, *args, **kwargs):
        use_cuda = kwargs.pop('use_cuda')
        n_samples = kwargs.pop('n_samples')
        latent, batch_indices, labels = get_latent(infer.model, data_loader, use_cuda=use_cuda)
        show_t_sne(latent, np.array([batch[0] for batch in batch_indices]), n_samples=n_samples)
