import time

import numpy as np

from scvi.metrics.classification import compute_accuracy
from scvi.metrics.log_likelihood import compute_log_likelihood
from scvi.models import VAE, VAEC, SVAEC


class Stats:
    def __init__(self, verbose=True, record_freq=5, n_epochs=-1, benchmark=False):
        self.benchmark = benchmark
        self.verbose = verbose
        self.record_freq = record_freq
        self.epoch = 0
        self.n_epochs = n_epochs
        self.begin = time.time()
        self.history = {"LL_train": [], "LL_test": [],
                        "Accuracy_train": [], "Accuracy_test": []}

    def callback(self, model, data_loader_train, data_loader_test, classifier=None):
        if self.epoch==0 or self.epoch==self.n_epochs or (self.epoch % self.record_freq == 0 and not self.benchmark):
            # In this case we do add the stats
            if self.verbose:
                print("EPOCH [%d/%d]: " % (self.epoch, self.n_epochs))
            self.add_ll_train(model, data_loader_train)
            self.add_ll_test(model, data_loader_test)
            self.add_accuracy_train(model, data_loader_train, classifier)
            self.add_accuracy_test(model, data_loader_test, classifier)
        self.epoch += 1

    def add_ll_train(self, model, data_loader):
        models = [VAE, VAEC, SVAEC]
        if type(model) in models:
            log_likelihood_train = compute_log_likelihood(model, data_loader)
            self.history["LL_train"].append(log_likelihood_train)
            if self.verbose:
                print("LL train is: %4f" % log_likelihood_train)

    def add_ll_test(self, model, data_loader):
        models = [VAE, VAEC, SVAEC]
        if type(model) in models:
            log_likelihood_test = compute_log_likelihood(model, data_loader)
            self.history["LL_test"].append(log_likelihood_test)
            if self.verbose:
                print("LL test is: %4f" % log_likelihood_test)

    def add_accuracy_train(self, model, data_loader, classifier=None):
        models = [VAEC, SVAEC]
        if type(model) in models or classifier:
            accuracy_train = compute_accuracy(model, data_loader, classifier)
            self.history["Accuracy_train"].append(accuracy_train)
            if self.verbose:
                print("Accuracy train is: %4f" % accuracy_train)

    def add_accuracy_test(self, model, data_loader, classifier=None):
        models = [VAEC, SVAEC]
        if type(model) in models or classifier:
            accuracy_test = compute_accuracy(model, data_loader, classifier)
            self.history["Accuracy_test"].append(accuracy_test)
            if self.verbose:
                print("Accuracy test is: %4f" % accuracy_test)

    def display_time(self):
        end = time.time()
        print("Total runtime for " + str(self.epoch) + " epochs is: " + str((end - self.begin))
              + " seconds for a mean per epoch runtime of " + str((end - self.begin) / self.epoch) + " seconds.")


class EarlyStopping:
    def __init__(self, patience=20, threshold=0.01, benchmark=False):
        self.benchmark = benchmark
        self.patience = patience
        self.threshold = threshold
        self.current_performances = np.ones((patience))

    def update(self, scalar):
        if self.benchmark:
            return True
        else:
            # Shift
            self.current_performances[:-1] = self.current_performances[1:]
            self.current_performances[-1] = scalar

            # Compute improvement
            improvement = ((self.current_performances[0] - self.current_performances[-1])
                           / self.current_performances[0])

            # Returns true if improvement is good enough
            return improvement >= self.threshold
