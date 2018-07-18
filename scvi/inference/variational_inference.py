import matplotlib.pyplot as plt
import numpy as np
import torch
from sklearn.manifold import TSNE
from torch.nn import functional as F

from scvi.dataset import CortexDataset
from scvi.dataset.data_loaders import DataLoaders
from scvi.dataset.data_loaders import TrainTestDataLoaders, AlternateSemiSupervisedDataLoaders, \
    JointSemiSupervisedDataLoaders
from scvi.metrics.classification import compute_accuracy, compute_accuracy_svc, compute_accuracy_rf
from scvi.metrics.clustering import get_latent, entropy_batch_mixing
from scvi.metrics.differential_expression import de_stats, de_cortex
from scvi.metrics.imputation import imputation
from scvi.metrics.log_likelihood import compute_log_likelihood
from . import Inference, ClassifierInference

plt.switch_backend('agg')


class VariationalInference(Inference):
    r"""The VariationalInference class for the unsupervised training of an autoencoder.

    Args:
        :model: A model instance from class ``VAE``, ``VAEC``, ``SVAEC``, ...
        :gene_dataset: A gene_dataset instance like ``CortexDataset()``
        :train_size: The train size, either a float between 0 and 1 or and integer for the number of training samples
         to use Default: ``0.8``.
        :**kwargs: Other keywords arguments from the general Inference class.

    Examples:
        >>> gene_dataset = CortexDataset()
        >>> vae = VAE(gene_dataset.nb_genes, n_batch=gene_dataset.n_batches * False,
        ... n_labels=gene_dataset.n_labels)

        >>> infer = VariationalInference(gene_dataset, vae, train_size=0.5)
        >>> infer.fit(n_epochs=20, lr=1e-3)
    """
    default_metrics_to_monitor = ['ll']

    def __init__(self, model, gene_dataset, train_size=0.8, **kwargs):
        super(VariationalInference, self).__init__(model, gene_dataset, **kwargs)
        self.kl = None
        self.data_loaders = TrainTestDataLoaders(self.gene_dataset, train_size=train_size, pin_memory=self.use_cuda)

    def loss(self, tensors):
        sample_batch, local_l_mean, local_l_var, batch_index, _ = tensors
        reconst_loss, kl_divergence = self.model(sample_batch, local_l_mean, local_l_var, batch_index)
        loss = torch.mean(reconst_loss + self.kl_weight * kl_divergence)
        return loss

    def on_epoch_begin(self):
        self.kl_weight = self.kl if self.kl is not None else min(1, self.epoch / self.n_epochs)

    def ll(self, name, verbose=False):
        ll = compute_log_likelihood(self.model, self.data_loaders[name], use_cuda=self.use_cuda)
        if verbose:
            print("LL for %s is : %.4f" % (name, ll))
        return ll

    ll.mode = 'min'

    def imputation_errors(self, name, **kwargs):
        return imputation(self.model, self.data_loaders[name], use_cuda=self.use_cuda, **kwargs)

    def imputation(self, name, verbose=False, *args, **kwargs):
        imputation_score = torch.median(self.imputation_errors(name, *args, **kwargs)).item()
        if verbose:
            print("Median Imputation Score for %s is : %.4f" % (name, imputation_score))
        return imputation_score

    imputation.mode = 'min'

    def differential_expression_stats(self, name, *args, **kwargs):
        return de_stats(self.model, self.data_loaders[name], *args, use_cuda=self.use_cuda, **kwargs)

    def differential_expression(self, name, *args, verbose=False, **kwargs):
        px_scale, all_labels = self.differential_expression_stats(name, *args, **kwargs)

        if type(self.gene_dataset) == CortexDataset:
            if 'use_cuda' in kwargs:
                kwargs.pop('use_cuda')
            de_score = de_cortex(px_scale, all_labels, self.gene_dataset.gene_names, **kwargs)
            if verbose:
                print("DE score for cortex on %s is : %.4f" % (name, de_score))
        return de_score

    differential_expression.mode = 'max'

    def entropy_batch_mixing(self, name, verbose=False, **kwargs):
        if self.gene_dataset.n_batches == 2:
            latent, batch_indices, labels = get_latent(self.model, self.data_loaders[name], use_cuda=self.use_cuda)
            be_score = entropy_batch_mixing(latent, batch_indices, **kwargs)
            if verbose:
                print("Entropy batch mixing :", be_score)
            return be_score

    entropy_batch_mixing.mode = 'max'

    def show_t_sne(self, name, n_samples=1000, color_by='', save_name=''):
        latent, batch_indices, labels = get_latent(self.model, self.data_loaders[name], use_cuda=self.use_cuda)
        idx_t_sne = np.random.permutation(len(latent))[:n_samples] if n_samples else np.arange(len(latent))
        if latent.shape[1] != 2:
            latent = TSNE().fit_transform(latent[idx_t_sne])
        if not color_by:
            plt.figure(figsize=(10, 10))
            plt.scatter(latent[:, 0], latent[:, 1])
        else:
            batch_indices = batch_indices[idx_t_sne].ravel()
            labels = labels[idx_t_sne].ravel()
            if color_by == 'batches' or color_by == 'labels':
                indices = batch_indices if color_by == 'batches' else labels
                n = self.gene_dataset.n_batches if color_by == 'batches' else self.gene_dataset.n_labels
                if hasattr(self.gene_dataset, 'cell_types') and color_by == 'labels':
                    plt_labels = self.gene_dataset.cell_types
                else:
                    plt_labels = [str(i) for i in range(len(np.unique(indices)))]
                plt.figure(figsize=(10, 10))
                for i, label in zip(range(n), plt_labels):
                    plt.scatter(latent[indices == i, 0], latent[indices == i, 1], label=label)
                plt.legend()
            elif color_by == 'batches and labels':
                fig, axes = plt.subplots(1, 2, figsize=(14, 7))
                for i in range(self.gene_dataset.n_batches):
                    axes[0].scatter(latent[batch_indices == i, 0], latent[batch_indices == i, 1], label=str(i))
                axes[0].set_title("batch coloring")
                axes[0].axis("off")
                axes[0].legend()

                indices = labels
                if hasattr(self.gene_dataset, 'cell_types'):
                    plt_labels = self.gene_dataset.cell_types
                else:
                    plt_labels = [str(i) for i in range(len(np.unique(indices)))]
                for i, cell_type in zip(range(self.gene_dataset.n_labels), plt_labels):
                    axes[1].scatter(latent[indices == i, 0], latent[indices == i, 1], label=cell_type)
                axes[1].set_title("label coloring")
                axes[1].axis("off")
                axes[1].legend()
        plt.axis("off")
        plt.tight_layout()
        if save_name:
            plt.savefig(save_name)


class SemiSupervisedVariationalInference(VariationalInference):
    r"""The abstract SemiSupervisedVariationalInference class for the semi-supervised training of an autoencoder.
    This parent class is inherited to specify the different training schemes for semi-supervised learning
    """
    default_metrics_to_monitor = VariationalInference.default_metrics_to_monitor + ['accuracy']

    def accuracy(self, name, verbose=False):
        acc = compute_accuracy(self.model, self.data_loaders[name], use_cuda=self.use_cuda)
        if verbose:
            print("Acc for %s is : %.4f" % (name, acc))
        return acc

    accuracy.mode = 'max'

    def svc_rf(self, **kwargs):
        if 'train' in self.data_loaders:
            raw_data = DataLoaders.raw_data(self.data_loaders['train'], self.data_loaders['test'])
        else:
            raw_data = DataLoaders.raw_data(self.data_loaders['labelled'], self.data_loaders['unlabelled'])
        (data_train, labels_train), (data_test, labels_test) = raw_data
        svc_scores = compute_accuracy_svc(data_train, labels_train, data_test, labels_test, **kwargs)
        rf_scores = compute_accuracy_rf(data_train, labels_train, data_test, labels_test, **kwargs)
        return svc_scores, rf_scores


class AlternateSemiSupervisedVariationalInference(SemiSupervisedVariationalInference):
    r"""The AlternateSemiSupervisedVariationalInference class for the semi-supervised training of an autoencoder.

    Args:
        :model: A model instance from class ``VAEC``, ``SVAEC``, ...
        :gene_dataset: A gene_dataset instance with pre-annotations like ``CortexDataset()``
        :n_labelled_samples_per_class: The number of labelled training samples per class. Default: ``50``.
        :**kwargs: Other keywords arguments from the general Inference class.

    Examples:
        >>> gene_dataset = CortexDataset()
        >>> svaec = SVAEC(gene_dataset.nb_genes, n_batch=gene_dataset.n_batches * False,
        ... n_labels=gene_dataset.n_labels)

        >>> infer = AlternateSemiSupervisedVariationalInference(gene_dataset, svaec, n_labelled_samples_per_class=10)
        >>> infer.fit(n_epochs=20, lr=1e-3)
    """
    def __init__(self, model, gene_dataset, n_labelled_samples_per_class=50, n_epochs_classifier=1,
                 lr_classification=0.1, **kwargs):
        super(AlternateSemiSupervisedVariationalInference, self).__init__(model, gene_dataset, **kwargs)

        self.n_epochs_classifier = n_epochs_classifier
        self.lr_classification = lr_classification
        self.data_loaders = AlternateSemiSupervisedDataLoaders(gene_dataset, n_labelled_samples_per_class)

        self.classifier_inference = ClassifierInference(
            model.classifier, gene_dataset, metrics_to_monitor=[], benchmark=True,
            data_loaders=self.data_loaders.classifier_data_loaders(), sampling_model=self.model
        )

    def on_epoch_end(self):
        self.classifier_inference.fit(self.n_epochs_classifier, lr=self.lr_classification)
        return super(AlternateSemiSupervisedVariationalInference, self).on_epoch_end()


class JointSemiSupervisedVariationalInference(SemiSupervisedVariationalInference):
    r"""The JointSemiSupervisedVariationalInference class for the semi-supervised training of an autoencoder.

    Args:
        :model: A model instance from class ``VAEC``, ``SVAEC``, ...
        :gene_dataset: A gene_dataset instance with pre-annotations like ``CortexDataset()``
        :n_labelled_samples_per_class: The number of labelled training samples per class. Default: ``50``.
        :**kwargs: Other keywords arguments from the general Inference class.

    Examples:
        >>> gene_dataset = CortexDataset()
        >>> svaec = SVAEC(gene_dataset.nb_genes, n_batch=gene_dataset.n_batches * False,
        ... n_labels=gene_dataset.n_labels)

        >>> infer = JointSemiSupervisedVariationalInference(gene_dataset, svaec, n_labelled_samples_per_class=10)
        >>> infer.fit(n_epochs=20, lr=1e-3)
    """
    def __init__(self, model, gene_dataset, n_labelled_samples_per_class=50, classification_ratio=100, **kwargs):
        super(JointSemiSupervisedVariationalInference, self).__init__(model, gene_dataset, **kwargs)
        self.data_loaders = JointSemiSupervisedDataLoaders(gene_dataset, n_labelled_samples_per_class)
        self.classification_ratio = classification_ratio

    def loss(self, tensors_all, tensors_labelled):
        loss = super(JointSemiSupervisedVariationalInference, self).loss(tensors_all)
        sample_batch, _, _, _, y = tensors_labelled
        classification_loss = F.cross_entropy(self.model.classify(sample_batch), y.view(-1))
        loss += classification_loss * self.classification_ratio
        return loss
