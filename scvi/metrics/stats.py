import time
from collections import defaultdict

import numpy as np

from scvi.metrics.classification import compute_accuracy
from scvi.metrics.log_likelihood import compute_log_likelihood
from scvi.models import VAE, VAEC, SVAEC


class Stats:
    def __init__(self, verbose=True, record_freq=5, n_epochs=-1, benchmark=False, names=['train', 'test', 'val']):
        self.verbose = verbose
        self.record_freq = record_freq
        self.n_epochs = n_epochs
        self.benchmark = benchmark
        self.names = names

        self.begin = time.time()
        self.epoch = 0
        self.history = defaultdict(lambda: [])

    def callback(self, model, *data_loaders, classifier=None):
        if self.epoch == 0 or self.epoch == self.n_epochs or (
                    self.epoch % self.record_freq == 0 and not self.benchmark):
            # In this case we do add the stats
            if self.verbose:
                print("EPOCH [%d/%d]: " % (self.epoch, self.n_epochs))
            for data_loader, name in zip(data_loaders, self.names):

                self.add_ll(model, data_loader, name=name)
                self.add_accuracy(model, data_loader, classifier=classifier, name=name)
        self.epoch += 1

    def add_ll(self, model, data_loader, name='train'):
        models = [VAE, VAEC, SVAEC]
        if type(model) in models:
            log_likelihood_test = compute_log_likelihood(model, data_loader)
            self.history["LL_%s" % name].append(log_likelihood_test)
            if self.verbose:
                print("LL %s is: %4f" % (name, log_likelihood_test))

    def add_accuracy(self, model, data_loader, classifier=None, name='train'):
        models = [VAEC, SVAEC]
        if type(model) in models or classifier:
            accuracy = compute_accuracy(model, data_loader, classifier)
            self.history["Accuracy_%s" % name].append(accuracy)
            if self.verbose:
                print("Accuracy %s is: %4f" % (name, accuracy))

    def display_time(self):
        end = time.time()
        print("Total runtime for " + str(self.epoch) + " epochs is: " + str((end - self.begin))
              + " seconds for a mean per epoch runtime of " + str((end - self.begin) / self.epoch) + " seconds.")


class EarlyStopping:
    def __init__(self, patience=250, threshold=0.01, benchmark=False):
        self.benchmark = benchmark
        self.patience = patience
        self.threshold = threshold
        self.current_performances = np.ones((patience))
        self.epoch = 0

    def update(self, scalar):
        self.epoch += 1
        if self.benchmark or self.epoch <= self.patience:
            return True
        else:
            # Shift
            self.current_performances[:-1] = self.current_performances[1:]
            self.current_performances[-1] = scalar

            # Compute improvement
            improvement = ((self.current_performances[-1] - self.current_performances[0])
                           / self.current_performances[0])

            # Returns true if improvement is good enough
            return improvement >= self.threshold
