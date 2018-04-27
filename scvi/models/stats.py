from scvi.utils import compute_accuracy
from scvi.log_likelihood import compute_log_likelihood
from . import VAE, VAEC, SVAEC


class Stats:
    def __init__(self, verbose=True, verbose_freq=5):
        self.verbose = verbose
        self.verbose_freq = verbose_freq
        self.n_epoch = 0
        self.history = {"LL_train": [], "LL_test": [],
                        "Accuracy_train": [], "Accuracy_test": []}

    def callback(self, model, data_loader_train, data_loader_test):
        if self.n_epoch % self.verbose_freq == 0:
            # In this case we do add the stats

            # Avoid dropout and batch normalization
            model.eval()

            print("For the epoch %d: " % self.n_epoch)
            self.add_ll_train(model, data_loader_train)
            self.add_ll_test(model, data_loader_test)
            self.add_accuracy_train(model, data_loader_train)
            self.add_accuracy_test(model, data_loader_test)
            model.train()
        self.n_epoch += 1

    def add_ll_train(self, model, data_loader):
        models = [VAE, VAEC, SVAEC]
        if type(model) in models:
            log_likelihood_train = compute_log_likelihood(model, data_loader)
            self.history["LL_train"].append(log_likelihood_train)
            print("LL train is: %4f" % log_likelihood_train)

    def add_ll_test(self, model, data_loader):
        models = [VAE, VAEC, SVAEC]
        if type(model) in models:
            log_likelihood_test = compute_log_likelihood(model, data_loader)
            self.history["LL_test"].append(log_likelihood_test)
            print("LL test is: %4f" % log_likelihood_test)

    def add_accuracy_train(self, model, data_loader):
        models = [VAEC, SVAEC]
        if type(model) in models:
            accuracy_train = compute_accuracy(model, data_loader)
            self.history["Accuracy_train"].append(accuracy_train)
            print("Accuracy train is: %4f" % accuracy_train)

    def add_accuracy_test(self, model, data_loader):
        models = [VAEC, SVAEC]
        if type(model) in models:
            accuracy_train = compute_accuracy(model, data_loader)
            self.history["Accuracy_test"].append(accuracy_train)
            print("Accuracy test is: %4f" % accuracy_train)
