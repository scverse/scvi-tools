import numpy as np
import torch

from scvi.dataset import CortexDataset
from scvi.dataset.utils import DataLoaderHandler
from scvi.metrics.adapt_encoder import adapt_encoder
from scvi.metrics.classification import compute_accuracy_rf, compute_accuracy_svc
from scvi.metrics.clustering import entropy_batch_mixing, get_latent
from scvi.metrics.differential_expression import de_stats, de_cortex
from scvi.metrics.imputation import imputation
from scvi.metrics.visualization import show_t_sne
from scvi.models import VAE
from scvi.train import Trainer, JointSemiSupervisedTrainer, AlternateSemiSupervisedTrainer


class Inference:
    def __init__(self, gene_dataset, model, use_cuda=True):
        self.gene_dataset = gene_dataset
        self.model = model
        self.use_cuda = use_cuda and torch.cuda.is_available()
        if self.use_cuda:
            self.model.cuda()

    def ll(self):
        if isinstance(self.model, VAE):
            best_ll = adapt_encoder(self.model, self.data_loader_test, n_path=1, n_epochs=1, record_freq=1,
                                    use_cuda=self.use_cuda)
            print("Best ll was :", best_ll)

        # - log-likelihood
        print("Log-likelihood Train:", self.stats.history["LL_train"][self.stats.best_index])
        print("Log-likelihood Test:", self.stats.history["LL_test"][self.stats.best_index])

    def imputation(self, rate=0.1):
        imputation_test = imputation(self.model, self.data_loader_test, rate=rate, use_cuda=self.use_cuda)
        print("Imputation score on test (MAE) is:", torch.median(imputation_test).item())

    def de(self):
        if type(self.gene_dataset) == CortexDataset:
            px_scale, all_labels = de_stats(self.model, self.data_loader_train, M_sampling=1, use_cuda=self.use_cuda)
            de_cortex(px_scale, all_labels, self.gene_dataset.gene_names, M_permutation=1)

    def show_tsne(self):
        if self.gene_dataset.n_batches == 2:
            latent, batch_indices, labels = get_latent(self.model, self.data_loader_train, use_cuda=self.use_cuda)
            print("Entropy batch mixing :", entropy_batch_mixing(latent, batch_indices))
            show_t_sne(latent, np.array([batch[0] for batch in batch_indices]))

    def all(self, unit_test=True):
        self.ll()
        self.imputation()
        self.de()
        self.show_tsne()


class UnsupervisedInference(Inference):
    def __init__(self, gene_dataset, model, use_cuda=True, train_size=0.1):
        super(UnsupervisedInference, self).__init__(gene_dataset, model, use_cuda=use_cuda)
        self.init_data_loaders(train_size)

    def init_data_loaders(self, train_size):
        self.data_loader_train, self.data_loader_test = (
            DataLoaderHandler(pin_memory=self.use_cuda).train_test_split(
                self.gene_dataset, train_size=train_size
            )
        )

    def train(self, n_epochs=1000, lr=1e-3, benchmark=False, **kargs):
        self.stats = Trainer().train(
            self.model, self.data_loader_train, self.data_loader_test, n_epochs=n_epochs, lr=lr, benchmark=benchmark,
            use_cuda=self.use_cuda, **kargs
        )
        return self.stats


class SemiSupervisedInference(Inference):
    def __init__(self, gene_dataset, model, use_cuda=True, n_labelled_samples_per_class=50):
        super(SemiSupervisedInference, self).__init__(gene_dataset, model, use_cuda=use_cuda)
        self.init_data_loaders(n_labelled_samples_per_class)

    def init_data_loaders(self, n_labelled_samples_per_class):
        self.data_loader_all, self.data_loader_labelled, self.data_loader_unlabelled = (
            DataLoaderHandler(pin_memory=self.use_cuda).all_labelled_unlabelled(
                self.gene_dataset, n_labelled_samples_per_class
            )
        )
        self.data_loader_train = self.data_loader_labelled
        self.data_loader_test = self.data_loader_unlabelled

    def train(self, n_epochs=1000, lr=1e-3, benchmark=False, mode="alternately", verbose=False, record_freq=5):
        if mode == "alternately":
            self.stats = AlternateSemiSupervisedTrainer().train(
                self.model, self.data_loader_all, self.data_loader_labelled, self.data_loader_unlabelled,
                n_epochs=n_epochs, lr=lr, benchmark=benchmark, use_cuda=self.use_cuda, verbose=verbose,
                record_freq=record_freq
            )
        elif mode == "jointly":
            self.stats = JointSemiSupervisedTrainer().train(
                self.model, self.data_loader_all, self.data_loader_labelled, self.data_loader_unlabelled,
                n_epochs=n_epochs, lr=lr, benchmark=benchmark, use_cuda=self.use_cuda, verbose=verbose,
                record_freq=record_freq
            )
        return self.stats

    def svc_rf(self, unit_test=False):
        (data_train, labels_train), (data_test, labels_test) = (
            DataLoaderHandler.raw_data(self.data_loader_labelled, self.data_loader_unlabelled)
        )
        accuracy_train_svc, accuracy_test_svc = compute_accuracy_svc(data_train, labels_train, data_test, labels_test,
                                                                     unit_test=unit_test)
        accuracy_train_rf, accuracy_test_rf = compute_accuracy_rf(data_train, labels_train, data_test, labels_test,
                                                                  unit_test=unit_test)
        return accuracy_train_svc, accuracy_test_svc, accuracy_train_rf, accuracy_test_rf

    def all(self, unit_test=False):
        super(SemiSupervisedInference, self).all()
        self.svc_rf(unit_test=unit_test)
