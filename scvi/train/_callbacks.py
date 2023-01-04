import warnings
from copy import deepcopy
from typing import Optional, Tuple

import numpy as np
import pytorch_lightning as pl
import torch
from pytorch_lightning.callbacks import Callback
from pytorch_lightning.callbacks.early_stopping import EarlyStopping
from pytorch_lightning.utilities import rank_zero_info

from scvi.dataloaders import AnnDataLoader


class SubSampleLabels(Callback):
    """Subsample labels."""

    def __init__(self):
        super().__init__()

    def on_train_epoch_start(self, trainer, pl_module):
        """Subsample labels at the beginning of each epoch."""
        trainer.train_dataloader.loaders.resample_labels()
        super().on_train_epoch_start(trainer, pl_module)


class SaveBestState(Callback):
    r"""
    Save the best module state and restore into model.

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

    def check_monitor_top(self, current):  # noqa: D102
        return self.monitor_op(current, self.best_module_metric_val)

    def on_val_epoch_end(self, trainer, pl_module):  # noqa: D102
        logs = trainer.callback_metrics
        self.epochs_since_last_check += 1

        if trainer.current_epoch > 0 and self.epochs_since_last_check >= self.period:
            self.epochs_since_last_check = 0
            current = logs.get(self.monitor)

            if current is None:
                warnings.warn(
                    f"Can save best module state only with {self.monitor} available,"
                    " skipping.",
                    RuntimeWarning,
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

    def on_train_start(self, trainer, pl_module):  # noqa: D102
        self.best_module_state = deepcopy(pl_module.module.state_dict())

    def on_train_end(self, trainer, pl_module):  # noqa: D102
        pl_module.module.load_state_dict(self.best_module_state)


class LoudEarlyStopping(EarlyStopping):
    """
    Wrapper of Pytorch Lightning EarlyStopping callback that prints the reason for stopping on teardown.

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

    def on_train_start(self, trainer, pl_module):  # noqa: D102
        module = pl_module.module
        if self.dataloader is None:
            dl = trainer.datamodule.train_dataloader()
        else:
            dl = self.dataloader
        module_init = module.init(module.rngs, next(iter(dl)))
        state, params = module_init.pop("params")
        pl_module.set_train_state(params, state)
