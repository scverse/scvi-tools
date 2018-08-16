import copy
import sys
import time
from collections import defaultdict
from itertools import cycle

import numpy as np
import torch
from sklearn.model_selection._split import _validate_shuffle_split
from torch.utils.data import DataLoader
from torch.utils.data.sampler import SequentialSampler, SubsetRandomSampler, RandomSampler
from tqdm import trange

from scvi.inference.posterior import Posterior

torch.set_grad_enabled(False)


class Trainer:
    r"""The abstract Inference class for training a PyTorch model and monitoring its statistics. It should be
    inherited at least with a .loss() function to be optimized in the training loop.

    Args:
        :model: A model instance from class ``VAE``, ``VAEC``, ``SVAEC``
        :gene_dataset: A gene_dataset instance like ``CortexDataset()``
        :use_cuda: Default: ``True``.
        :metrics_to_monitor: A list of the metrics to monitor. If not specified, will use the
            ``default_metrics_to_monitor`` as specified in each . Default: ``None``.
        :benchmark: if True, prevents statistics computation in the training. Default: ``False``.
        :verbose: If statistics should be displayed along training. Default: ``None``.
        :frequency: The frequency at which to keep track of statistics. Default: ``None``.
        :early_stopping_metric: The statistics on which to perform early stopping. Default: ``None``.
        :save_best_state_metric:  The statistics on which we keep the network weights achieving the best store, and
            restore them at the end of training. Default: ``None``.
        :on: The data_loader name reference for the ``early_stopping_metric`` and ``save_best_state_metric``, that
            should be specified if any of them is. Default: ``None``.
    """
    default_metrics_to_monitor = []

    def __init__(self, model, gene_dataset, use_cuda=True, metrics_to_monitor=None, benchmark=False,
                 verbose=False, frequency=None, early_stopping_metric=None,
                 save_best_state_metric=None, on=None, early_stopping_kwargs={'patience': 15, 'threshold': 3},
                 weight_decay=1e-6, data_loaders_kwargs=dict()):

        self.model = model
        self.gene_dataset = gene_dataset

        self.data_loaders_kwargs = {
            "batch_size": 128,
            "pin_memory": use_cuda
        }
        self.data_loaders_kwargs.update(data_loaders_kwargs)

        self.weight_decay = weight_decay
        self.benchmark = benchmark
        self.epoch = 0
        self.training_time = 0

        if metrics_to_monitor is not None:
            self.metrics_to_monitor = metrics_to_monitor
        else:
            self.metrics_to_monitor = self.default_metrics_to_monitor

        self.early_stopping_metric = early_stopping_metric
        self.save_best_state_metric = save_best_state_metric
        self.on = on
        mode = getattr(Posterior, early_stopping_metric).mode if early_stopping_metric is not None else None
        mode_save_state = getattr(Posterior,
                                  save_best_state_metric).mode if save_best_state_metric is not None else None
        self.early_stopping = EarlyStopping(benchmark=benchmark, mode=mode, mode_save_state=mode_save_state,
                                            **early_stopping_kwargs)

        self.use_cuda = use_cuda and torch.cuda.is_available()
        if self.use_cuda:
            self.model.cuda()

        self.frequency = frequency if not benchmark else None
        self.verbose = verbose

        self.history = defaultdict(lambda: [])
        self.posteriors_dict = dict()
        self.posteriors_loop = []

    def data_loaders_loop(self):  # returns an zipped iterable corresponding to loss signature
        data_loaders_loop = [self.posteriors_dict[name].data_loader for name in self.posteriors_loop]
        return zip(data_loaders_loop[0], *[cycle(data_loader) for data_loader in data_loaders_loop[1:]])

    def compute_metrics(self):
        begin = time.time()
        with torch.set_grad_enabled(False):
            self.model.eval()
            if self.frequency and (
                            self.epoch == 0 or self.epoch == self.n_epochs or (self.epoch % self.frequency == 0)):
                if self.verbose:
                    print("\nEPOCH [%d/%d]: " % (self.epoch, self.n_epochs))

                for name, posterior in self.posteriors_dict.items():
                    if hasattr(posterior, 'to_monitor'):
                        for metric in posterior.to_monitor:
                            result = getattr(posterior, metric)(verbose=self.verbose)
                            self.history[metric + '_' + name] += [result]
            self.model.train()
            self.compute_metrics_time += time.time() - begin

    def train(self, n_epochs=20, lr=1e-3, params=None):
        begin = time.time()
        with torch.set_grad_enabled(True):
            self.model.train()

            if params is None:
                params = filter(lambda p: p.requires_grad, self.model.parameters())

            optimizer = torch.optim.Adam(params, lr=lr, weight_decay=self.weight_decay)
            self.epoch = 0
            self.compute_metrics_time = 0
            self.n_epochs = n_epochs
            self.compute_metrics()

            with trange(n_epochs, desc="training", file=sys.stdout, disable=self.verbose) as pbar:
                # We have to use tqdm this way so it works in Jupyter notebook.
                # See https://stackoverflow.com/questions/42212810/tqdm-in-jupyter-notebook
                self.on_epoch_begin()

                for epoch in pbar:
                    pbar.update(1)
                    for tensors_list in self.data_loaders_loop():
                        loss = self.loss(*tensors_list)
                        optimizer.zero_grad()
                        loss.backward()
                        optimizer.step()

                    if not self.on_epoch_end():
                        break

            if self.save_best_state_metric is not None:
                self.model.load_state_dict(self.best_state_dict)
                self.compute_metrics()

            self.model.eval()
            self.training_time += (time.time() - begin) - self.compute_metrics_time
            if self.verbose and self.frequency:
                print("\nTraining time:  %i s. / %i epochs" % (int(self.training_time), self.n_epochs))

    def on_epoch_begin(self):
        pass

    def on_epoch_end(self):
        self.epoch += 1
        self.compute_metrics()
        if self.save_best_state_metric is not None and self.on is not None:
            if self.early_stopping.update_state(self.history[self.save_best_state_metric + '_' + self.on][-1]):
                self.best_state_dict = self.model.state_dict()
                self.best_epoch = self.epoch

        continue_training = True
        if self.early_stopping_metric is not None and self.on is not None:
            continue_training = self.early_stopping.update(
                self.history[self.early_stopping_metric + '_' + self.on][-1]
            )
        return continue_training

    def train_test(self, model=None, gene_dataset=None, train_size=0.1, test_size=None, seed=0, suffix='',
                   cls=Posterior):
        """
        :param train_size: float, int, or None (default is 0.1)
        :param test_size: float, int, or None (default is None)
        """
        model = self.model if model is None and hasattr(self, "model") else model
        gene_dataset = self.gene_dataset if gene_dataset is None and hasattr(self, "model") else gene_dataset
        n = len(gene_dataset)
        n_train, n_test = _validate_shuffle_split(n, test_size, train_size)
        np.random.seed(seed=seed)
        permutation = np.random.permutation(n)
        indices_test = permutation[:n_test]
        indices_train = permutation[n_test:(n_test + n_train)]

        return (
            self.create_posterior(model, gene_dataset, indices=indices_train, name='train' + suffix, cls=cls),
            self.create_posterior(model, gene_dataset, indices=indices_test, name='test' + suffix, cls=cls)
        )

    def create_posterior(self, model=None, gene_dataset=None, shuffle=False, indices=None, name='', cls=Posterior):
        model = self.model if model is None and hasattr(self, "model") else model
        gene_dataset = self.gene_dataset if gene_dataset is None and hasattr(self, "model") else gene_dataset
        self.data_loaders_kwargs.update({'collate_fn': gene_dataset.collate_fn})
        if indices is not None and shuffle:
            raise ValueError('indices is mutually exclusive with shuffle')
        if indices is None:
            if shuffle:
                sampler = RandomSampler(gene_dataset)
            else:
                sampler = SequentialSampler(gene_dataset)
        else:
            sampler = SubsetRandomSampler(indices)
        return cls(model, gene_dataset,
                   ToCudaDataLoader(gene_dataset, use_cuda=self.use_cuda, sampler=sampler, **self.data_loaders_kwargs),
                   name=name)


class SequentialSubsetSampler(SubsetRandomSampler):
    def __init__(self, indices):
        self.indices = np.sort(indices)

    def __iter__(self):
        return iter(self.indices)


class ToCudaDataLoader(DataLoader):
    def __init__(self, dataset, use_cuda=True, **kwargs):
        self.kwargs = kwargs
        self.use_cuda = use_cuda and torch.cuda.is_available()
        super(ToCudaDataLoader, self).__init__(dataset, **kwargs)

    def to_cuda(self, tensors):
        return [t.cuda(async=self.use_cuda) if self.use_cuda else t for t in tensors]

    def sequential(self, batch_size=128):
        kwargs = copy.copy(self.kwargs)
        kwargs['batch_size'] = batch_size
        if hasattr(self, 'sampler') and hasattr(self.sampler, 'indices'):
            kwargs['sampler'] = SequentialSubsetSampler(indices=self.sampler.indices)
        else:
            kwargs['sampler'] = SequentialSampler(self.dataset)
        return ToCudaDataLoader(self.dataset, use_cuda=self.use_cuda, **kwargs)

    def __iter__(self):
        return map(self.to_cuda, super(ToCudaDataLoader, self).__iter__())


class EarlyStopping:
    def __init__(self, patience=15, threshold=3, benchmark=False, mode="max", mode_save_state="max"):
        self.benchmark = benchmark
        self.patience = patience
        self.threshold = threshold
        self.epoch = 0
        self.wait = 0
        self.mode = mode
        # We set the best to + inf because we're dealing with a loss we want to minimize
        self.current_performance = np.inf
        self.best_performance = np.inf
        self.best_performance_state = np.inf
        # If we want to maximize, we start at - inf
        if self.mode == "max":
            self.best_performance *= -1
            self.current_performance *= -1
        self.mode_save_state = mode_save_state
        if self.mode_save_state == "max":
            self.best_performance_state *= -1

    def update(self, scalar):
        self.epoch += 1
        if self.benchmark or self.epoch < self.patience:
            continue_training = True
        elif self.wait >= self.patience:
            continue_training = False
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

            continue_training = True
        if not continue_training:
            print("\nStopping early: no improvement of more than " + str(self.threshold) +
                  " nats in " + str(self.patience) + " epochs")
            print("If the early stopping criterion is too strong, "
                  "please instantiate it with different parameters in the train method.")
        return continue_training

    def update_state(self, scalar):
        improved = ((self.mode_save_state == "max" and scalar - self.best_performance_state > 0) or
                    (self.mode_save_state == "min" and self.best_performance_state - scalar > 0))
        if improved:
            self.best_performance_state = scalar
        return improved
