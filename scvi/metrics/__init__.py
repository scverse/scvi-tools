import numpy as np
import torch

from scvi.dataset import CortexDataset
from scvi.dataset.utils import DataLoaders
from scvi.metrics.clustering import get_latent, entropy_batch_mixing
from .classification import compute_accuracy, compute_accuracy_svc, compute_accuracy_rf
from .differential_expression import de_stats, de_cortex
from .imputation import imputation
from .log_likelihood import compute_log_likelihood
from .visualization import show_t_sne


class VariationalInferencePlugin:
    def ll(self, name):
        ll = compute_log_likelihood(self.model, self.data_loaders[name], use_cuda=self.use_cuda)
        print("LL for %s is : %.4f" % (name, ll))
        return ll

    ll.mode = 'min'

    def imputation_stats(self, name, **kwargs):
        return imputation(self.model, self.data_loaders[name], use_cuda=self.use_cuda, **kwargs)

    def imputation(self, infer, data_loader, name, *args, **kwargs):
        imputation_score = torch.median(self.imputation_task(infer, data_loader, name, *args, **kwargs)).item()
        print("Median Imputation Score for %s is : %.4f" % (name, imputation_score))
        return imputation_score

    imputation.mode = 'min'

    def differential_expression_stats(self, name, *args, **kwargs):
        return de_stats(self.model, self.data_loaders[name], *args, use_cuda=self.use_cuda, **kwargs)

    def differential_expression(self, name, *args, **kwargs):
        px_scale, all_labels = self.differential_expression_stats(name, *args, **kwargs)

        if type(self.gene_dataset) == CortexDataset:
            if 'use_cuda' in kwargs:
                kwargs.pop('use_cuda')
            de_score = de_cortex(px_scale, all_labels, self.gene_dataset.gene_names, **kwargs)
            print("DE score for cortex on %s is : %.4f" % (name, de_score))
        return de_score

    differential_expression.mode = 'max'

    def batch_entropy_mixing(self, name, *args, **kwargs):
        if self.gene_dataset.n_batches == 2:
            latent, batch_indices, labels = get_latent(self.model, self.data_loaders[name], use_cuda=self.use_cuda)
            be_score = entropy_batch_mixing(latent, batch_indices, **kwargs)
            print("Entropy batch mixing :", be_score)
            return be_score

    batch_entropy_mixing.mode = 'max'

    def show_t_sne(self, name, **kwargs):
        n_samples = kwargs.pop('n_samples')
        latent, batch_indices, labels = get_latent(self.model, self.data_loaders[name], use_cuda=self.use_cuda)
        show_t_sne(latent, np.array([batch[0] for batch in batch_indices]), n_samples=n_samples)


class AccuracyPlugin:
    def accuracy(self, name, *args, **kwargs):
        acc = compute_accuracy(self, self.data_loaders[name], *args, use_cuda=self.use_cuda, **kwargs)
        print("Acc for %s is : %.4f" % (name, acc))
        return acc

    accuracy.mode = 'max'

    def svc_rf(self, **kwargs):
        if 'train' in self.data_loaders:
            raw_data = DataLoaders.raw_data(self.data_loaders['train'], self.data_loaders['test'])
        else:
            raw_data = DataLoaders.raw_data(self.data_loaders['labelled'], self.data_loaders['unlabelled'])
        (data_train, labels_train), (data_test, labels_test) = raw_data
        score_train = compute_accuracy_svc(data_train, labels_train, data_test, labels_test, **kwargs)
        score_test = compute_accuracy_rf(data_train, labels_train, data_test, labels_test, **kwargs)
        return score_train, score_test


class SemiSupervisedVariationalInferencePlugin(VariationalInferencePlugin, AccuracyPlugin):
    pass
