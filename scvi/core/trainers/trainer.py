import logging
import time
from abc import abstractmethod
from collections import OrderedDict, defaultdict
from itertools import cycle
from typing import List

import anndata
import numpy as np
import torch
from sklearn.model_selection._split import _validate_shuffle_split
from torch.utils.data.sampler import SubsetRandomSampler

from scvi import _CONSTANTS
from scvi._utils import track
from scvi.core.data_loaders import ScviDataLoader

logger = logging.getLogger(__name__)


class Trainer:
    """
    The abstract Trainer class for training a PyTorch model and monitoring its statistics.

    It should be inherited at least with a ``.loss()`` function to be optimized in the training loop.

    Parameters
    ----------
    model :
        A model instance from class ``VAE``, ``VAEC``, ``SCANVI``
    adata:
        A registered anndata object
    use_cuda :
        Default: ``True``.
    metrics_to_monitor :
        A list of the metrics to monitor. If not specified, will use the
        ``default_metrics_to_monitor`` as specified in each . Default: ``None``.
    benchmark :
        if True, prevents statistics computation in the training. Default: ``False``.
    frequency :
        The frequency at which to keep track of statistics. Default: ``None``.
    early_stopping_metric :
        The statistics on which to perform early stopping. Default: ``None``.
    save_best_state_metric :
        The statistics on which we keep the network weights achieving the best store, and
        restore them at the end of training. Default: ``None``.
    on :
        The data_loader name reference for the ``early_stopping_metric`` and ``save_best_state_metric``, that
        should be specified if any of them is. Default: ``None``.
    silent :
        If True, disables progress bar.
    seed :
        Random seed for train/test/validate split
    """

    default_metrics_to_monitor = []

    def __init__(
        self,
        model,
        adata: anndata.AnnData,
        use_cuda: bool = True,
        metrics_to_monitor: List = None,
        benchmark: bool = False,
        frequency: int = None,
        weight_decay: float = 1e-6,
        early_stopping_kwargs: dict = None,
        data_loader_kwargs: dict = None,
        silent: bool = False,
        batch_size: int = 128,
        seed: int = 0,
        max_nans: int = 10,
    ):

        # Model, dataset management
        self.model = model
        self.adata = adata
        self._scvi_data_loaders = OrderedDict()
        self.seed = seed  # For train/test splitting
        self.use_cuda = use_cuda and torch.cuda.is_available()
        if self.use_cuda:
            self.model.cuda()

        # Data loader attributes
        self.batch_size = batch_size
        self.data_loader_kwargs = {"pin_memory": use_cuda}
        data_loader_kwargs = data_loader_kwargs if data_loader_kwargs else dict()
        self.data_loader_kwargs.update(data_loader_kwargs)

        # Optimization attributes
        self.optimizer = None
        self.weight_decay = weight_decay
        self.n_epochs = None
        self.epoch = -1  # epoch = self.epoch + 1 in compute metrics
        self.training_time = 0
        self.n_iter = 0

        # Training NaNs handling
        self.max_nans = max_nans
        self.current_loss = None  # torch.Tensor training loss
        self.previous_loss_was_nan = False
        self.nan_counter = 0  # Counts occuring NaNs during training

        # Metrics and early stopping
        self.compute_metrics_time = None
        if metrics_to_monitor is not None:
            self.metrics_to_monitor = set(metrics_to_monitor)
        else:
            self.metrics_to_monitor = set(self.default_metrics_to_monitor)
        early_stopping_kwargs = (
            early_stopping_kwargs if early_stopping_kwargs else dict()
        )
        self.early_stopping = EarlyStopping(**early_stopping_kwargs)
        self.benchmark = benchmark
        self.frequency = frequency if not benchmark else None
        self.history = defaultdict(list)
        self.best_state_dict = self.model.state_dict()
        self.best_epoch = self.epoch

        self.silent = silent

    @torch.no_grad()
    def compute_metrics(self):
        computed = set()  # for ensuring no double computation in function
        begin = time.time()
        epoch = self.epoch + 1
        if self.frequency and (
            epoch == 0 or epoch == self.n_epochs or (epoch % self.frequency == 0)
        ):
            with torch.set_grad_enabled(False):
                self.model.eval()
                logger.debug("\nEPOCH [%d/%d]: " % (epoch, self.n_epochs))

                for name, scdl in self._scvi_data_loaders.items():
                    message = " ".join([s.capitalize() for s in name.split("_")[-2:]])
                    if scdl.n_cells < 5:
                        logging.debug(
                            message + " is too small to track metrics (<5 samples)"
                        )
                        continue
                    if hasattr(scdl, "to_monitor"):
                        for metric in scdl.to_monitor:
                            if metric not in self.metrics_to_monitor:
                                logger.debug(message)
                                result = getattr(scdl, metric)()
                                out_str = metric + "_" + name
                                self.history[out_str] += [result]
                                computed.add(out_str)
                    for metric in self.metrics_to_monitor:
                        result = getattr(scdl, metric)()
                        out_str = metric + "_" + name
                        self.history[out_str] += [result]
                        computed.add(out_str)
                self.model.train()
        # compute metrics every epoch if using early stopping
        if self.early_stopping.early_stopping_metric is not None:
            name = self.early_stopping.on
            metric = self.early_stopping.early_stopping_metric
            out_str = metric + "_" + name
            if out_str not in computed:
                scdl = self._scvi_data_loaders[name]
                result = getattr(scdl, metric)()
                self.history[metric + "_" + name] += [result]

        self.compute_metrics_time += time.time() - begin

    def train(self, n_epochs=400, lr=1e-3, eps=0.01, params=None, **extras_kwargs):
        begin = time.time()
        self.model.train()

        if params is None:
            params = filter(lambda p: p.requires_grad, self.model.parameters())

        self.optimizer = torch.optim.Adam(
            params, lr=lr, eps=eps, weight_decay=self.weight_decay
        )

        # Initialization of other model's optimizers
        self.training_extras_init(**extras_kwargs)

        self.compute_metrics_time = 0
        self.n_epochs = n_epochs
        self.compute_metrics()

        self.on_training_begin()

        for self.epoch in track(
            range(n_epochs), description="Training...", disable=self.silent
        ):
            self.on_epoch_begin()
            for tensors_dict in self.data_loaders_loop():
                if tensors_dict[0][_CONSTANTS.X_KEY].shape[0] < 3:
                    continue
                self.on_iteration_begin()
                # Update the model's parameters after seeing the data
                self.on_training_loop(tensors_dict)
                # Checks the training status, ensures no nan loss
                self.on_iteration_end()

            # Computes metrics and controls early stopping
            if not self.on_epoch_end():
                break
        if self.early_stopping.save_best_state_metric is not None:
            self.model.load_state_dict(self.best_state_dict)
            self.compute_metrics()

        self.model.eval()
        self.training_extras_end()

        self.training_time += (time.time() - begin) - self.compute_metrics_time
        self.on_training_end()

    def on_training_loop(self, tensors_dict):
        self.current_loss = loss = self.loss(*tensors_dict)
        self.optimizer.zero_grad()
        loss.backward()
        self.optimizer.step()

    def training_extras_init(self, **extras_kwargs):
        """Other necessary models to simultaneously train."""
        pass

    def training_extras_end(self):
        """Place to put extra models in eval mode, etc."""
        pass

    def on_training_begin(self):
        pass

    def on_epoch_begin(self):
        # Epochs refer to a pass through the entire dataset (in minibatches)
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
                logger.info("Reducing LR on epoch {}.".format(self.epoch))
                for param_group in self.optimizer.param_groups:
                    param_group["lr"] *= self.early_stopping.lr_factor

        return continue_training

    def on_iteration_begin(self):
        # Iterations refer to minibatches
        pass

    def on_iteration_end(self):
        self.check_training_status()
        self.n_iter += 1

    def on_training_end(self):
        pass

    def check_training_status(self):
        """
        Checks if loss is admissible.

        If not, training is stopped after max_nans consecutive inadmissible loss
        loss corresponds to the training loss of the model.

        `max_nans` is the maximum number of consecutive NaNs after which a ValueError will be
        """
        loss_is_nan = torch.isnan(self.current_loss).item()
        if loss_is_nan:
            logger.warning("Model training loss was NaN")
            self.nan_counter += 1
            self.previous_loss_was_nan = True
        else:
            self.nan_counter = 0
            self.previous_loss_was_nan = False

        if self.nan_counter >= self.max_nans:
            raise ValueError(
                "Loss was NaN {} consecutive times: the model is not training properly. "
                "Consider using a lower learning rate.".format(self.max_nans)
            )

    @property
    @abstractmethod
    def scvi_data_loaders_loop(self):
        pass

    def data_loaders_loop(self):
        """Returns an zipped iterable corresponding to loss signature."""
        data_loaders_loop = [
            self._scvi_data_loaders[name] for name in self.scvi_data_loaders_loop
        ]
        return zip(
            data_loaders_loop[0],
            *[cycle(data_loader) for data_loader in data_loaders_loop[1:]]
        )

    def register_data_loader(self, name, value):
        name = name.strip("_")
        self._scvi_data_loaders[name] = value

    def __getattr__(self, name):
        if "_scvi_data_loaders" in self.__dict__:
            _scvi_data_loaders = self.__dict__["_scvi_data_loaders"]
            if name.strip("_") in _scvi_data_loaders:
                return _scvi_data_loaders[name.strip("_")]
        return object.__getattribute__(self, name)

    def __delattr__(self, name):
        if name.strip("_") in self._scvi_data_loaders:
            del self._scvi_data_loaders[name.strip("_")]
        else:
            object.__delattr__(self, name)

    def __setattr__(self, name, value):
        if isinstance(value, ScviDataLoader):
            name = name.strip("_")
            self.register_data_loader(name, value)
        else:
            object.__setattr__(self, name, value)

    def train_test_validation(
        self,
        model=None,
        adata=None,
        train_size=0.9,
        test_size=None,
        type_class=ScviDataLoader,
    ):
        """
        Creates data loaders ``train_set``, ``test_set``, ``validation_set``.

        If ``train_size + test_size < 1`` then ``validation_set`` is non-empty.

        Parameters
        ----------
        train_size :
            float, or None (default is 0.9)
        test_size :
            float, or None (default is None)
        model :
             (Default value = None)
        adata:
             (Default value = None)
        type_class :
             (Default value = ScviDataLoader)
        """
        train_size = float(train_size)
        if train_size > 1.0 or train_size <= 0.0:
            raise ValueError(
                "train_size needs to be greater than 0 and less than or equal to 1"
            )

        model = self.model if model is None and hasattr(self, "model") else model
        # what to do here if it has the attribute modle
        adata = self.adata if adata is None and hasattr(self, "model") else adata
        n = len(adata)
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
            self.create_scvi_dl(
                model, adata, indices=indices_train, type_class=type_class
            ),
            self.create_scvi_dl(
                model, adata, indices=indices_test, type_class=type_class
            ),
            self.create_scvi_dl(
                model, adata, indices=indices_validation, type_class=type_class
            ),
        )

    def create_scvi_dl(
        self,
        model=None,
        adata: anndata.AnnData = None,
        shuffle=False,
        indices=None,
        type_class=ScviDataLoader,
    ):
        model = self.model if model is None and hasattr(self, "model") else model
        adata = self.adata if adata is None and hasattr(self, "model") else adata
        return type_class(
            model,
            adata,
            shuffle=shuffle,
            indices=indices,
            use_cuda=self.use_cuda,
            batch_size=self.batch_size,
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
        scvi_data_loader_class=ScviDataLoader,
    ):
        self.benchmark = benchmark
        self.patience = patience
        self.threshold = threshold
        self.epoch = 0
        self.wait = 0
        self.wait_lr = 0
        self.mode = (
            getattr(scvi_data_loader_class, early_stopping_metric).mode
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
            getattr(ScviDataLoader, save_best_state_metric).mode
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
