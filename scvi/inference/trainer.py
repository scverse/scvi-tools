import logging
import sys
import time

from abc import abstractmethod
from collections import defaultdict, OrderedDict
from itertools import cycle

import numpy as np
import torch

from sklearn.model_selection._split import _validate_shuffle_split
from torch.utils.data.sampler import SubsetRandomSampler
from tqdm import trange

from scvi.inference.posterior import Posterior

logger = logging.getLogger(__name__)


class Trainer:
    r"""The abstract Trainer class for training a PyTorch model and monitoring its statistics. It should be
    inherited at least with a .loss() function to be optimized in the training loop.

    Args:
        :model: A model instance from class ``VAE``, ``VAEC``, ``SCANVI``
        :gene_dataset: A gene_dataset instance like ``CortexDataset()``
        :use_cuda: Default: ``True``.
        :metrics_to_monitor: A list of the metrics to monitor. If not specified, will use the
            ``default_metrics_to_monitor`` as specified in each . Default: ``None``.
        :benchmark: if True, prevents statistics computation in the training. Default: ``False``.
        :frequency: The frequency at which to keep track of statistics. Default: ``None``.
        :early_stopping_metric: The statistics on which to perform early stopping. Default: ``None``.
        :save_best_state_metric:  The statistics on which we keep the network weights achieving the best store, and
            restore them at the end of training. Default: ``None``.
        :on: The data_loader name reference for the ``early_stopping_metric`` and ``save_best_state_metric``, that
            should be specified if any of them is. Default: ``None``.
        :show_progbar: If False, disables progress bar.
        :seed: Random seed for train/test/validate split
    """
    default_metrics_to_monitor = []

    def __init__(
        self,
        model,
        gene_dataset,
        use_cuda=True,
        metrics_to_monitor=None,
        benchmark=False,
        frequency=None,
        weight_decay=1e-6,
        early_stopping_kwargs=None,
        data_loader_kwargs=None,
        show_progbar=True,
        seed=0,
    ):
        # handle mutable defaults
        early_stopping_kwargs = (
            early_stopping_kwargs if early_stopping_kwargs else dict()
        )
        data_loader_kwargs = data_loader_kwargs if data_loader_kwargs else dict()

        self.model = model
        self.gene_dataset = gene_dataset
        self._posteriors = OrderedDict()
        self.seed = seed

        self.data_loader_kwargs = {"batch_size": 128, "pin_memory": use_cuda}
        self.data_loader_kwargs.update(data_loader_kwargs)

        self.weight_decay = weight_decay
        self.benchmark = benchmark
        self.epoch = -1  # epoch = self.epoch + 1 in compute metrics
        self.training_time = 0

        if metrics_to_monitor is not None:
            self.metrics_to_monitor = set(metrics_to_monitor)
        else:
            self.metrics_to_monitor = set(self.default_metrics_to_monitor)

        self.early_stopping = EarlyStopping(**early_stopping_kwargs)

        if self.early_stopping.early_stopping_metric:
            self.metrics_to_monitor.add(self.early_stopping.early_stopping_metric)

        self.use_cuda = use_cuda and torch.cuda.is_available()
        if self.use_cuda:
            self.model.cuda()

        self.frequency = frequency if not benchmark else None

        self.history = defaultdict(list)

        self.best_state_dict = self.model.state_dict()
        self.best_epoch = self.epoch

        self.show_progbar = show_progbar

    @torch.no_grad()
    def compute_metrics(self):
        begin = time.time()
        epoch = self.epoch + 1
        if self.frequency and (
            epoch == 0 or epoch == self.n_epochs or (epoch % self.frequency == 0)
        ):
            with torch.set_grad_enabled(False):
                self.model.eval()
                logger.debug("\nEPOCH [%d/%d]: " % (epoch, self.n_epochs))

                for name, posterior in self._posteriors.items():
                    message = " ".join([s.capitalize() for s in name.split("_")[-2:]])
                    if posterior.nb_cells < 5:
                        logging.debug(
                            message + " is too small to track metrics (<5 samples)"
                        )
                        continue
                    if hasattr(posterior, "to_monitor"):
                        for metric in posterior.to_monitor:
                            if metric not in self.metrics_to_monitor:
                                logger.debug(message)
                                result = getattr(posterior, metric)()
                                self.history[metric + "_" + name] += [result]
                    for metric in self.metrics_to_monitor:
                        result = getattr(posterior, metric)()
                        self.history[metric + "_" + name] += [result]
                self.model.train()
        self.compute_metrics_time += time.time() - begin

    def train(self, n_epochs=20, lr=1e-3, eps=0.01, params=None):
        begin = time.time()
        self.model.train()

        if params is None:
            params = filter(lambda p: p.requires_grad, self.model.parameters())

        optimizer = self.optimizer = torch.optim.Adam(
            params, lr=lr, eps=eps, weight_decay=self.weight_decay
        )

        self.compute_metrics_time = 0
        self.n_epochs = n_epochs
        self.compute_metrics()

        with trange(
            n_epochs, desc="training", file=sys.stdout, disable=not self.show_progbar
        ) as pbar:
            # We have to use tqdm this way so it works in Jupyter notebook.
            # See https://stackoverflow.com/questions/42212810/tqdm-in-jupyter-notebook
            for self.epoch in pbar:
                self.on_epoch_begin()
                pbar.update(1)
                for tensors_list in self.data_loaders_loop():
                    if tensors_list[0][0].shape[0] < 3:
                        continue
                    loss = self.loss(*tensors_list)
                    optimizer.zero_grad()
                    loss.backward()
                    optimizer.step()

                if not self.on_epoch_end():
                    break

        if self.early_stopping.save_best_state_metric is not None:
            self.model.load_state_dict(self.best_state_dict)
            self.compute_metrics()

        self.model.eval()
        self.training_time += (time.time() - begin) - self.compute_metrics_time
        if self.frequency:
            logger.debug(
                "\nTraining time:  %i s. / %i epochs"
                % (int(self.training_time), self.n_epochs)
            )

    def on_epoch_begin(self):
        pass

    def on_epoch_end(self):
        self.compute_metrics()
        on = self.early_stopping.on
        early_stopping_metric = self.early_stopping.early_stopping_metric
        save_best_state_metric = self.early_stopping.save_best_state_metric
        if save_best_state_metric is not None and on is not None:
            if self.early_stopping.update_state(
                self.history[save_best_state_metric + "_" + on][-1]
            ):
                self.best_state_dict = self.model.state_dict()
                self.best_epoch = self.epoch

        continue_training = True
        if early_stopping_metric is not None and on is not None:
            continue_training, reduce_lr = self.early_stopping.update(
                self.history[early_stopping_metric + "_" + on][-1]
            )
            if reduce_lr:
                logger.info("Reducing LR.")
                for param_group in self.optimizer.param_groups:
                    param_group["lr"] *= self.early_stopping.lr_factor

        return continue_training

    @property
    @abstractmethod
    def posteriors_loop(self):
        pass

    def data_loaders_loop(self):
        """returns an zipped iterable corresponding to loss signature"""
        data_loaders_loop = [self._posteriors[name] for name in self.posteriors_loop]
        return zip(
            data_loaders_loop[0],
            *[cycle(data_loader) for data_loader in data_loaders_loop[1:]]
        )

    def register_posterior(self, name, value):
        name = name.strip("_")
        self._posteriors[name] = value

    def corrupt_posteriors(
        self, rate=0.1, corruption="uniform", update_corruption=True
    ):
        if not hasattr(self.gene_dataset, "corrupted") and update_corruption:
            self.gene_dataset.corrupt(rate=rate, corruption=corruption)
        for name, posterior in self._posteriors.items():
            self.register_posterior(name, posterior.corrupted())

    def uncorrupt_posteriors(self):
        for name_, posterior in self._posteriors.items():
            self.register_posterior(name_, posterior.uncorrupted())

    def __getattr__(self, name):
        if "_posteriors" in self.__dict__:
            _posteriors = self.__dict__["_posteriors"]
            if name.strip("_") in _posteriors:
                return _posteriors[name.strip("_")]
        return object.__getattribute__(self, name)

    def __delattr__(self, name):
        if name.strip("_") in self._posteriors:
            del self._posteriors[name.strip("_")]
        else:
            object.__delattr__(self, name)

    def __setattr__(self, name, value):
        if isinstance(value, Posterior):
            name = name.strip("_")
            self.register_posterior(name, value)
        else:
            object.__setattr__(self, name, value)

    def train_test_validation(
        self,
        model=None,
        gene_dataset=None,
        train_size=0.1,
        test_size=None,
        type_class=Posterior,
    ):
        """Creates posteriors ``train_set``, ``test_set``, ``validation_set``.
            If ``train_size + test_size < 1`` then ``validation_set`` is non-empty.

            :param train_size: float, int, or None (default is 0.1)
            :param test_size: float, int, or None (default is None)
            """
        model = self.model if model is None and hasattr(self, "model") else model
        gene_dataset = (
            self.gene_dataset
            if gene_dataset is None and hasattr(self, "model")
            else gene_dataset
        )
        n = len(gene_dataset)
        try:
            n_train, n_test = _validate_shuffle_split(n, test_size, train_size)
        except ValueError:
            if train_size != 1.0:
                raise ValueError(
                    "Choice of train_size={} and test_size={} not understood".format(
                        train_size, test_size
                    )
                )
            n_train, n_test = n, 0
        random_state = np.random.RandomState(seed=self.seed)
        permutation = random_state.permutation(n)
        indices_test = permutation[:n_test]
        indices_train = permutation[n_test : (n_test + n_train)]
        indices_validation = permutation[(n_test + n_train) :]

        return (
            self.create_posterior(
                model, gene_dataset, indices=indices_train, type_class=type_class
            ),
            self.create_posterior(
                model, gene_dataset, indices=indices_test, type_class=type_class
            ),
            self.create_posterior(
                model, gene_dataset, indices=indices_validation, type_class=type_class
            ),
        )

    def create_posterior(
        self,
        model=None,
        gene_dataset=None,
        shuffle=False,
        indices=None,
        type_class=Posterior,
    ):
        model = self.model if model is None and hasattr(self, "model") else model
        gene_dataset = (
            self.gene_dataset
            if gene_dataset is None and hasattr(self, "model")
            else gene_dataset
        )
        return type_class(
            model,
            gene_dataset,
            shuffle=shuffle,
            indices=indices,
            use_cuda=self.use_cuda,
            data_loader_kwargs=self.data_loader_kwargs,
        )


class SequentialSubsetSampler(SubsetRandomSampler):
    def __init__(self, indices):
        self.indices = np.sort(indices)

    def __iter__(self):
        return iter(self.indices)


class EarlyStopping:
    def __init__(
        self,
        early_stopping_metric: str = None,
        save_best_state_metric: str = None,
        on: str = "test_set",
        patience: int = 15,
        threshold: int = 3,
        benchmark: bool = False,
        reduce_lr_on_plateau: bool = False,
        lr_patience: int = 10,
        lr_factor: float = 0.5,
        posterior_class=Posterior,
    ):
        self.benchmark = benchmark
        self.patience = patience
        self.threshold = threshold
        self.epoch = 0
        self.wait = 0
        self.wait_lr = 0
        self.mode = (
            getattr(posterior_class, early_stopping_metric).mode
            if early_stopping_metric is not None
            else None
        )
        # We set the best to + inf because we're dealing with a loss we want to minimize
        self.current_performance = np.inf
        self.best_performance = np.inf
        self.best_performance_state = np.inf
        # If we want to maximize, we start at - inf
        if self.mode == "max":
            self.best_performance *= -1
            self.current_performance *= -1
        self.mode_save_state = (
            getattr(Posterior, save_best_state_metric).mode
            if save_best_state_metric is not None
            else None
        )
        if self.mode_save_state == "max":
            self.best_performance_state *= -1

        self.early_stopping_metric = early_stopping_metric
        self.save_best_state_metric = save_best_state_metric
        self.on = on
        self.reduce_lr_on_plateau = reduce_lr_on_plateau
        self.lr_patience = lr_patience
        self.lr_factor = lr_factor

    def update(self, scalar):
        self.epoch += 1
        if self.benchmark:
            continue_training = True
            reduce_lr = False
        elif self.wait >= self.patience:
            continue_training = False
            reduce_lr = False
        else:
            # Check if we should reduce the learning rate
            if not self.reduce_lr_on_plateau:
                reduce_lr = False
            elif self.wait_lr >= self.lr_patience:
                reduce_lr = True
                self.wait_lr = 0
            else:
                reduce_lr = False
            # Shift
            self.current_performance = scalar

            # Compute improvement
            if self.mode == "max":
                improvement = self.current_performance - self.best_performance
            elif self.mode == "min":
                improvement = self.best_performance - self.current_performance
            else:
                raise NotImplementedError("Unknown optimization mode")

            # updating best performance
            if improvement > 0:
                self.best_performance = self.current_performance

            if improvement < self.threshold:
                self.wait += 1
                self.wait_lr += 1
            else:
                self.wait = 0
                self.wait_lr = 0

            continue_training = True
        if not continue_training:
            # FIXME: log total number of epochs run
            logger.info(
                "\nStopping early: no improvement of more than "
                + str(self.threshold)
                + " nats in "
                + str(self.patience)
                + " epochs"
            )
            logger.info(
                "If the early stopping criterion is too strong, "
                "please instantiate it with different parameters in the train method."
            )
        return continue_training, reduce_lr

    def update_state(self, scalar):
        improved = (
            self.mode_save_state == "max" and scalar - self.best_performance_state > 0
        ) or (
            self.mode_save_state == "min" and self.best_performance_state - scalar > 0
        )
        if improved:
            self.best_performance_state = scalar
        return improved
