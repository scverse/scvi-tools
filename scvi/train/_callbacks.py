import warnings
from copy import deepcopy
from typing import Callable, Dict, List, Optional, Tuple, Union

import flax
import lightning.pytorch as pl
import numpy as np
import torch
from lightning.pytorch.callbacks import Callback
from lightning.pytorch.callbacks.early_stopping import EarlyStopping
from lightning.pytorch.utilities import rank_zero_info

from scvi import settings
from scvi.dataloaders import AnnDataLoader
from scvi.model.base import BaseModelClass

MetricCallable = Callable[[BaseModelClass], float]


class MetricsCallback(Callback):
    """Computes metrics on validation end and logs them to the logger.

    Parameters
    ----------
    metric_fns
        Validation metrics to compute and log. One of the following:

        * `:class:`~scvi.train._callbacks.MetricCallable`: A function that takes in a
            :class:`~scvi.model.base.BaseModelClass` and returns a `float`.
            The function's `__name__`is used as the logging name.

        * `List[:class:`~scvi.train._callbacks.MetricCallable`]`: Same as above but in
            a list.

        * `Dict[str, :class:`~scvi.train._callbacks.MetricCallable`]`: Same as above,
            but the keys are used as the logging names instead.
    """

    def __init__(
        self,
        metric_fns: Union[
            MetricCallable, List[MetricCallable], Dict[str, MetricCallable]
        ],
    ):
        super().__init__()

        if callable(metric_fns):
            metric_fns = [metric_fns]

        if not isinstance(metric_fns, (list, dict)):
            raise TypeError("`metric_fns` must be a `list` or `dict`.")

        values = metric_fns if isinstance(metric_fns, list) else metric_fns.values()
        for val in values:
            if not callable(val):
                raise TypeError("`metric_fns` must contain functions only.")

        if not isinstance(metric_fns, dict):
            metric_fns = {f.__name__: f for f in metric_fns}

        self.metric_fns = metric_fns

    def on_validation_end(self, trainer, pl_module):
        """Compute metrics at the end of validation.

        Sets the model to trained mode before computing metrics and restores training
        mode thereafter. Metrics are not logged with a `"validation"` prefix as the
        metrics are only computed on the validation set.
        """
        model = trainer._model  # TODO: Remove with a better way to access model
        model.is_trained = True

        metrics = {}
        for metric_name, metric_fn in self.metric_fns.items():
            metric_value = metric_fn(model)
            metrics[metric_name] = metric_value

        metrics["epoch"] = trainer.current_epoch

        pl_module.logger.log_metrics(metrics, trainer.global_step)
        model.is_trained = False


class SubSampleLabels(Callback):
    """Subsample labels."""

    def __init__(self):
        super().__init__()

    def on_train_epoch_start(self, trainer, pl_module):
        """Subsample labels at the beginning of each epoch."""
        trainer.train_dataloader.resample_labels()
        super().on_train_epoch_start(trainer, pl_module)


class SaveBestState(Callback):
    r"""Save the best module state and restore into model.

    Parameters
    ----------
    monitor
        quantity to monitor.
    verbose
        verbosity, True or False.
    mode
        one of ["min", "max"].
    period
        Interval (number of epochs) between checkpoints.

    Examples
    --------
    from scvi.train import Trainer
    from scvi.train import SaveBestState
    """

    def __init__(
        self,
        monitor: str = "elbo_validation",
        mode: str = "min",
        verbose=False,
        period=1,
    ):
        super().__init__()

        self.monitor = monitor
        self.verbose = verbose
        self.period = period
        self.epochs_since_last_check = 0
        self.best_module_state = None

        if mode not in ["min", "max"]:
            raise ValueError(
                f"SaveBestState mode {mode} is unknown",
            )

        if mode == "min":
            self.monitor_op = np.less
            self.best_module_metric_val = np.Inf
            self.mode = "min"
        elif mode == "max":
            self.monitor_op = np.greater
            self.best_module_metric_val = -np.Inf
            self.mode = "max"
        else:
            if "acc" in self.monitor or self.monitor.startswith("fmeasure"):
                self.monitor_op = np.greater
                self.best_module_metric_val = -np.Inf
                self.mode = "max"
            else:
                self.monitor_op = np.less
                self.best_module_metric_val = np.Inf
                self.mode = "min"

    def check_monitor_top(self, current):
        return self.monitor_op(current, self.best_module_metric_val)

    def on_validation_epoch_end(self, trainer, pl_module):
        logs = trainer.callback_metrics
        self.epochs_since_last_check += 1

        if trainer.current_epoch > 0 and self.epochs_since_last_check >= self.period:
            self.epochs_since_last_check = 0
            current = logs.get(self.monitor)

            if current is None:
                warnings.warn(
                    f"Can save best module state only with {self.monitor} available, "
                    "skipping.",
                    RuntimeWarning,
                    stacklevel=settings.warnings_stacklevel,
                )
            else:
                if isinstance(current, torch.Tensor):
                    current = current.item()
                if self.check_monitor_top(current):
                    self.best_module_state = deepcopy(pl_module.module.state_dict())
                    self.best_module_metric_val = current

                    if self.verbose:
                        rank_zero_info(
                            f"\nEpoch {trainer.current_epoch:05d}: {self.monitor} reached."
                            f" Module best state updated."
                        )

    def on_train_start(self, trainer, pl_module):
        self.best_module_state = deepcopy(pl_module.module.state_dict())

    def on_train_end(self, trainer, pl_module):
        pl_module.module.load_state_dict(self.best_module_state)


class LoudEarlyStopping(EarlyStopping):
    """Wrapper of Pytorch Lightning EarlyStopping callback that prints the reason for stopping on teardown.

    When the early stopping condition is met, the reason is saved to the callback instance,
    then printed on teardown. By printing on teardown, we do not interfere with the progress
    bar callback.
    """

    def __init__(self, **kwargs) -> None:
        super().__init__(**kwargs)
        self.early_stopping_reason = None

    def _evaluate_stopping_criteria(self, current: torch.Tensor) -> Tuple[bool, str]:
        should_stop, reason = super()._evaluate_stopping_criteria(current)
        if should_stop:
            self.early_stopping_reason = reason
        return should_stop, reason

    def teardown(
        self,
        _trainer: pl.Trainer,
        _pl_module: pl.LightningModule,
        stage: Optional[str] = None,
    ) -> None:
        """Print the reason for stopping on teardown."""
        if self.early_stopping_reason is not None:
            print(self.early_stopping_reason)


class JaxModuleInit(Callback):
    """A callback to initialize the Jax-based module."""

    def __init__(self, dataloader: AnnDataLoader = None) -> None:
        super().__init__()
        self.dataloader = dataloader

    def on_train_start(self, trainer, pl_module):
        module = pl_module.module
        if self.dataloader is None:
            dl = trainer.datamodule.train_dataloader()
        else:
            dl = self.dataloader
        module_init = module.init(module.rngs, next(iter(dl)))
        state, params = flax.core.pop(module_init, "params")
        pl_module.set_train_state(params, state)
