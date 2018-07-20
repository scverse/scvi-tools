import sys
from collections import defaultdict

import torch
from tqdm import trange

from scvi.metrics.early_stopping import EarlyStopping

torch.set_grad_enabled(False)


class Inference:
    r"""The abstract Inference class for training a PyTorch model and monitoring its statistics. It should be
    inherited at least with a .loss() function to be optimized in the training loop.

    Args:
        :model: A model instance from class ``VAE``, ``VAEC``, ``SVAEC``, ...
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

    def __init__(self, model, gene_dataset, use_cuda=True, metrics_to_monitor=None, data_loaders=None, benchmark=False,
                 verbose=False, frequency=None, early_stopping_metric=None, save_best_state_metric=None, on=None):
        self.model = model
        self.gene_dataset = gene_dataset
        self.data_loaders = data_loaders
        self.benchmark = benchmark
        self.epoch = 0

        if metrics_to_monitor is not None:
            self.metrics_to_monitor = metrics_to_monitor
        else:
            self.metrics_to_monitor = self.default_metrics_to_monitor

        self.early_stopping_metric = early_stopping_metric
        self.save_best_state_metric = save_best_state_metric
        self.on = on
        mode = getattr(self, early_stopping_metric).mode if early_stopping_metric is not None else None
        mode_save_state = getattr(self, save_best_state_metric).mode if save_best_state_metric is not None else None
        self.early_stopping = EarlyStopping(benchmark=benchmark, mode=mode, mode_save_state=mode_save_state)

        self.use_cuda = use_cuda and torch.cuda.is_available()
        if self.use_cuda:
            self.model.cuda()

        self.frequency = frequency if not benchmark else None
        self.verbose = verbose

        self.history = defaultdict(lambda: [])

    def compute_metrics(self):
        with torch.set_grad_enabled(False):
            self.model.eval()
            if self.frequency and (
                        self.epoch == 0 or self.epoch == self.n_epochs or (self.epoch % self.frequency == 0)):
                if self.verbose:
                    print("\nEPOCH [%d/%d]: " % (self.epoch, self.n_epochs))
                for name in self.data_loaders.to_monitor:
                    for metric in self.metrics_to_monitor:
                        result = getattr(self, metric)(name=name, verbose=self.verbose)
                        self.history[metric + '_' + name] += [result]
            self.model.train()

    def fit(self, n_epochs=20, lr=1e-3):
        with torch.set_grad_enabled(True):
            self.model.train()
            optimizer = self.optimizer if hasattr(self, 'optimizer') else \
                torch.optim.Adam(filter(lambda p: p.requires_grad, self.model.parameters()), lr=lr)
            self.epoch = 0
            self.n_epochs = n_epochs
            self.compute_metrics()
            with trange(n_epochs, desc="training", file=sys.stdout, disable=self.verbose) as pbar:
                # We have to use tqdm this way so it works in Jupyter notebook.
                # See https://stackoverflow.com/questions/42212810/tqdm-in-jupyter-notebook
                self.on_epoch_begin()
                for epoch in pbar:
                    pbar.update(1)
                    for tensors_list in self.data_loaders:
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
