import sys
from collections import defaultdict

import torch
# model, *data_loaders are attributes of this class
from tqdm import trange

from scvi.metrics import AccuracyMetric, LLMetric
# from scvi.inference.inferencetask import EarlyStopping
from scvi.metrics.early_stopping import EarlyStopping
from scvi.utils import to_cuda, enable_grad

metrics = {'ll': LLMetric(),
           'accuracy': AccuracyMetric()
           }


class Inference:
    """
    Inference class cares about the general training of a PyTorch model with its statistics
    """
    metrics = []

    def __init__(self, model, gene_dataset, data_loaders=None, metrics_to_monitor=[], use_cuda=True, benchmark=False,
                 record_freq=1, verbose=True, early_stopping_metric=None, save_best_state_metric=None, on=None):
        self.model = model
        self.gene_dataset = gene_dataset
        self.data_loaders = data_loaders
        self.benchmark = benchmark
        self.epoch = 0

        self.metrics_to_monitor = metrics_to_monitor
        self.early_stopping_metric = early_stopping_metric
        self.save_best_state_metric = save_best_state_metric
        self.on = on
        mode = metrics[early_stopping_metric].mode if early_stopping_metric is not None else None
        mode_save_state = metrics[save_best_state_metric].mode if save_best_state_metric is not None else None
        self.early_stopping = EarlyStopping(benchmark=benchmark, mode=mode, mode_save_state=mode_save_state)
        self.use_cuda = use_cuda and torch.cuda.is_available()
        if self.use_cuda:
            self.model.cuda()
        self.set_metrics_methods()

        self.verbose = verbose if not benchmark else False
        self.record_freq = record_freq

        self.history = defaultdict(lambda: [])

    def set_metrics_methods(self):
        def metric_method_getter(metric):
            def metric_method(name, **kargs):
                return metrics[metric](self, self.data_loaders[name], name=name, use_cuda=self.use_cuda, **kargs)

            return metric_method

        for metric in self.metrics:
            setattr(self, metric, metric_method_getter(metric))

    def compute_metrics(self):
        if self.verbose and (self.epoch == 0 or self.epoch == self.n_epochs or (self.epoch % self.record_freq == 0)):
            print("EPOCH [%d/%d]: " % (self.epoch, self.n_epochs))
            for name in self.data_loaders.to_monitor:
                for metric in self.metrics_to_monitor:
                    result = getattr(self, metric)(name=name)
                    self.history[metric + '_' + name] += [result]

    @enable_grad()
    def fit(self, n_epochs=20, lr=1e-3):
        optimizer = torch.optim.Adam(filter(lambda p: p.requires_grad, self.model.parameters()), lr=lr)
        self.epoch = 0
        self.n_epochs = n_epochs
        self.compute_metrics()
        # Training the model

        with trange(n_epochs, desc="training", file=sys.stdout, disable=self.benchmark or self.verbose) as pbar:
            # We have to use tqdm this way so it works in Jupyter notebook.
            # See https://stackoverflow.com/questions/42212810/tqdm-in-jupyter-notebook
            self.on_epoch_begin()
            for epoch in pbar:
                pbar.update(1)
                for tensors_list in self.data_loaders:
                    loss = self.loss(*[to_cuda(tensors, use_cuda=self.use_cuda) for tensors in tensors_list])
                    optimizer.zero_grad()
                    loss.backward()
                    optimizer.step()
                if not self.on_epoch_end():
                    break
        if self.save_best_state_metric is not None:
            self.model.load_state_dict(self.best_state_dict)
            print("Now re-computing metric")
            self.compute_metrics()

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
            continue_training = self.early_stopping.update(self.history[self.early_stopping_metric + '_' + self.on][-1])
        return continue_training
