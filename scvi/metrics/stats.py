import time
from collections import defaultdict

import numpy as np

from scvi.metrics.classification import compute_accuracy
from scvi.metrics.log_likelihood import compute_log_likelihood
from scvi.models import VAE, VAEC, SVAEC


class Stats:
    def __init__(self, verbose=True, record_freq=5, n_epochs=-1, benchmark=False, names=['train', 'test', 'val'],
                 use_cuda=True):
        self.verbose = verbose if not benchmark else False
        self.record_freq = record_freq
        self.n_epochs = n_epochs
        self.benchmark = benchmark
        self.names = names

        self.best_params = None
        self.best_score = +np.inf
        self.begin = time.time()
        self.epoch = 0
        self.history = defaultdict(lambda: [])
        self.use_cuda = use_cuda

    def callback(self, model, *data_loaders, classifier=None):
        if self.epoch == 0 or self.epoch == self.n_epochs or (
                        self.epoch % self.record_freq == 0 and not self.benchmark):
            # In this case we do add the stats
            if self.verbose:
                print("EPOCH [%d/%d]: " % (self.epoch, self.n_epochs))
            for data_loader, name in zip(data_loaders, self.names):
                self.add_ll(model, data_loader, name=name)
                self.add_accuracy(model, data_loader, classifier=classifier, name=name)
            self.save_best_params(model)
        self.epoch += 1

    def add_ll(self, model, data_loader, name='train'):
        available_models = [VAE, VAEC, SVAEC]
        if type(model) in available_models:
            log_likelihood = compute_log_likelihood(model, data_loader, use_cuda=self.use_cuda)
            self.history["LL_%s" % name].append(log_likelihood)
            if self.verbose:
                print("LL %s is: %4f" % (name, log_likelihood))

    def add_accuracy(self, model, data_loader, classifier=None, name='train'):
        available_models = [VAEC, SVAEC]
        if type(model) in available_models or classifier:
            accuracy = compute_accuracy(model, data_loader, classifier, use_cuda=self.use_cuda)
            self.history["Accuracy_%s" % name].append(accuracy)
            if self.verbose:
                print("Accuracy %s is: %4f" % (name, accuracy))

    def save_best_params(self, model):
        current_score = self.history["LL_%s" % self.names[0]][-1]
        if current_score < self.best_score:
            self.best_index = self.epoch // self.record_freq
            self.best_score = current_score
            self.best_params = model.state_dict()

    def set_best_params(self, model):
        model.load_state_dict(self.best_params)

    def display_time(self):
        end = time.time()
        print("Total runtime for " + str(self.epoch) + " epochs is: " + str((end - self.begin))
              + " seconds for a mean per epoch runtime of " + str((end - self.begin) / self.epoch) + " seconds.")


class EarlyStopping:
    def __init__(self, patience=15, threshold=3, benchmark=False, mode="min"):
        self.benchmark = benchmark
        self.patience = patience
        self.threshold = threshold
        self.epoch = 0
        self.wait = 0
        self.mode = mode
        # We set the best to + inf because we're dealing with a loss we want to minimize
        self.current_performance = np.inf
        self.best_performance = np.inf
        # If we want to maximize, we start at - inf
        if self.mode == "max":
            self.best *= -1
            self.current_performance *= -1

    def update(self, scalar):
        self.epoch += 1
        if self.benchmark or self.epoch < self.patience:
            return True
        if self.wait >= self.patience:
            return False
        else:
            # Shift
            self.current_performance = scalar

            # Compute improvement
            if self.mode == "max":
                improvement = self.current_performance - self.best_performance
            elif self.mode == "min":
                improvement = self.best_performance - self.current_performance

            # updating best performance
            if improvement > 0:
                self.best_performance = self.current_performance

            if improvement < self.threshold:
                self.wait += 1
            else:
                self.wait = 0

            return True
