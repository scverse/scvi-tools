from __future__ import annotations

import logging
import os
import warnings
from collections.abc import Callable
from copy import deepcopy
from datetime import datetime
from shutil import rmtree
from typing import TYPE_CHECKING

import flax
import numpy as np
import pyro
import torch
from lightning.pytorch.callbacks import Callback, ModelCheckpoint
from lightning.pytorch.callbacks.early_stopping import EarlyStopping
from lightning.pytorch.utilities import rank_zero_info
from lightning.pytorch.utilities.rank_zero import rank_prefixed_message

from scvi import settings
from scvi.model.base import BaseModelClass
from scvi.model.base._save_load import _load_saved_files

if TYPE_CHECKING:
    import lightning.pytorch as pl

    from scvi.dataloaders import AnnDataLoader

MetricCallable = Callable[[BaseModelClass], float]

log = logging.getLogger(__name__)


class SaveCheckpoint(ModelCheckpoint):
    """``BETA`` Saves model checkpoints based on a monitored metric.

    Inherits from :class:`~lightning.pytorch.callbacks.ModelCheckpoint` and modifies the default
    behavior to save the full model state instead of just the state dict. This is necessary for
    compatibility with :class:`~scvi.model.base.BaseModelClass`.

    The best model save and best model score based on ``monitor`` can be accessed post-training
    with the ``best_model_path`` and ``best_model_score`` attributes, respectively.

    Known issues:

    * Does not set ``train_indices``, ``validation_indices``, and ``test_indices`` for checkpoints.
    * Does not set ``history`` for checkpoints. This can be accessed in the final model however.
    * Unsupported arguments: ``save_weights_only`` and ``save_last``.

    Parameters
    ----------
    dirpath
        Base directory to save the model checkpoints. If ``None``, defaults to a subdirectory in
        :attr:``scvi.settings.logging_dir`` formatted with the current date, time, and monitor.
    filename
        Name for the checkpoint directories, which can contain formatting options for auto-filling.
        If ``None``, defaults to ``{epoch}-{step}-{monitor}``.
    monitor
        Metric to monitor for checkpointing.
    load_best_on_end
        If ``True``, loads the best model state into the model at the end of training.
    check_nan_gradients
        If ``True``, will use the on exception callback to store best model in case of training
        exception caused by NaN's in gradients or loss calculations.
    **kwargs
        Additional keyword arguments passed into the constructor for
        :class:`~lightning.pytorch.callbacks.ModelCheckpoint`.
    """

    def __init__(
        self,
        dirpath: str | None = None,
        filename: str | None = None,
        monitor: str = "validation_loss",
        load_best_on_end: bool = False,
        check_nan_gradients: bool = False,
        **kwargs,
    ):
        if dirpath is None:
            dirpath = os.path.join(
                settings.logging_dir,
                datetime.now().strftime("%Y-%m-%d_%H-%M-%S") + f"_{monitor}",
            )
        if filename is None:
            filename = "{epoch}-{step}-{" + monitor + "}"
        if "save_weights_only" in kwargs:
            warnings.warn(
                "`save_weights_only` is not supported in `SaveCheckpoint` and will be ignored.",
                RuntimeWarning,
                stacklevel=settings.warnings_stacklevel,
            )
            kwargs.pop("save_weights_only")
        if "save_last" in kwargs:
            warnings.warn(
                "`save_last` is not supported in `SaveCheckpoint` and will be ignored.",
                RuntimeWarning,
                stacklevel=settings.warnings_stacklevel,
            )
            kwargs.pop("save_last")
        self.load_best_on_end = load_best_on_end
        self.loss_is_nan = False
        self.check_nan_gradients = check_nan_gradients

        super().__init__(
            dirpath=dirpath,
            filename=filename,
            monitor=monitor,
            **kwargs,
        )

    def on_save_checkpoint(self, trainer: pl.Trainer, *args) -> None:
        """Saves the model state on Lightning checkpoint saves."""
        # set post training state before saving
        model = trainer._model
        model.module.eval()
        model.is_trained_ = True
        model.trainer = trainer

        monitor_candidates = self._monitor_candidates(trainer)
        save_path = self.format_checkpoint_name(monitor_candidates)
        # by default, the function above gives a .ckpt extension
        save_path = save_path.split(".ckpt")[0]
        model.save(save_path, save_anndata=False, overwrite=True)

        model.module.train()
        model.is_trained_ = False

    def _remove_checkpoint(self, trainer: pl.Trainer, filepath: str) -> None:
        """Removes model saves that are no longer needed.

        Calls the superclass method and then removes the :class:`~scvi.model.base.BaseModelClass`
        save directory.
        """
        super()._remove_checkpoint(trainer, filepath)

        model_path = filepath.split(".ckpt")[0]
        if os.path.exists(model_path) and os.path.isdir(model_path):
            rmtree(model_path)

    def _update_best_and_save(
        self,
        current: torch.Tensor,
        trainer: pl.Trainer,
        monitor_candidates: dict[str, torch.Tensor],
    ) -> None:
        """Replaces Lightning checkpoints with :class:`~scvi.model.base.BaseModelClass` saves.

        Calls the superclass method and then replaces the Lightning checkpoint file with
        the :class:`~scvi.model.base.BaseModelClass` save directory.
        """
        super()._update_best_and_save(current, trainer, monitor_candidates)

        if os.path.exists(self.best_model_path):
            os.remove(self.best_model_path)
        self.best_model_path = self.best_model_path.split(".ckpt")[0]

    def on_train_end(self, trainer: pl.Trainer, pl_module: pl.LightningModule) -> None:
        """Loads the best model state into the model at the end of training."""
        if not self.load_best_on_end:
            return

        _, _, best_state_dict, _ = _load_saved_files(
            self.best_model_path,
            load_adata=False,
            map_location=pl_module.module.device,
        )
        pyro_param_store = best_state_dict.pop("pyro_param_store", None)
        pl_module.module.load_state_dict(best_state_dict)
        if pyro_param_store is not None:
            # For scArches shapes are changed and we don't want to overwrite these changed shapes.
            pyro.get_param_store().set_state(pyro_param_store)

    def on_train_batch_end(self, trainer, pl_module, outputs, batch, batch_idx) -> None:
        # Check for NaN in the loss
        loss = outputs.get("loss") if isinstance(outputs, dict) else outputs
        if torch.isnan(loss).any():
            self.loss_is_nan = True

    def on_exception(self, trainer, pl_module, exception) -> None:
        """Save the model in case of unexpected exceptions, like Nan in loss or gradients"""
        if (not isinstance(exception, KeyboardInterrupt)) and self.check_nan_gradients:
            if self.loss_is_nan:
                self.reason = (
                    "\033[31m[Warning] NaN detected in the loss. Stopping training. "
                    "Please verify your model and data. "
                    "Saving model....Please load it back and continue training\033[0m"
                )
            else:
                self.reason = (
                    "\033[31m[Warning] Exception occurred during training (Nan or Inf gradients). "
                    "Please verify your model and data. "
                    "Saving model....Please load it back and continue training\033[0m"
                )
            trainer.should_stop = True

            # Loads the best model state into the model after the exception.
            _, _, best_state_dict, _ = _load_saved_files(
                self.best_model_path,
                load_adata=False,
                map_location=pl_module.module.device,
            )
            pyro_param_store = best_state_dict.pop("pyro_param_store", None)
            pl_module.module.load_state_dict(best_state_dict)
            if pyro_param_store is not None:
                # For scArches shapes are changed and we don't want to overwrite
                # these changed shapes.
                pyro.get_param_store().set_state(pyro_param_store)
            print(self.reason + f". Model saved to {self.best_model_path}")
            self._log_info(trainer, self.reason, False)
            return

    @staticmethod
    def _log_info(trainer: pl.Trainer, message: str, log_rank_zero_only: bool) -> None:
        rank = trainer.global_rank if trainer.world_size > 1 else None
        message = rank_prefixed_message(message, rank)
        if rank is None or not log_rank_zero_only or rank == 0:
            log.info(message)


class SubSampleLabels(Callback):
    """Subsample labels."""

    def __init__(self):
        super().__init__()

    def on_train_epoch_start(self, trainer, pl_module):
        """Subsample labels at the beginning of each epoch."""
        trainer.train_dataloader.resample_labels()
        super().on_train_epoch_start(trainer, pl_module)


class SaveBestState(Callback):
    r"""``DEPRECATED`` Save the best module state and restore into model.

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

    Notes
    -----
    Lifecycle: deprecated in v1.2 and to be removed in v1.3. Please use
        :class:`~scvi.train.callbacks.SaveCheckpoint` instead.
    """

    def __init__(
        self,
        monitor: str = "elbo_validation",
        mode: str = "min",
        verbose=False,
        period=1,
    ):
        super().__init__()

        warnings.warn(
            "`SaveBestState` is deprecated in v1.2 and will be removed in v1.3. Please use "
            "`SaveCheckpoint` instead. See https://github.com/scverse/scvi-tools/issues/2568 "
            "for more details.",
            DeprecationWarning,
            stacklevel=settings.warnings_stacklevel,
        )

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
            self.best_module_metric_val = np.inf
            self.mode = "min"
        elif mode == "max":
            self.monitor_op = np.greater
            self.best_module_metric_val = -np.inf
            self.mode = "max"
        else:
            if "acc" in self.monitor or self.monitor.startswith("fmeasure"):
                self.monitor_op = np.greater
                self.best_module_metric_val = -np.inf
                self.mode = "max"
            else:
                self.monitor_op = np.less
                self.best_module_metric_val = np.inf
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
                    f"Can save best module state only with {self.monitor} available, skipping.",
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
    """Loud early stopping callback.

    Wrapper of :class:`~lightning.pytorch.callbacks.early_stopping.EarlyStopping callback that
    prints the reason for stopping on teardown. When the early stopping condition is met, the
    reason is saved to the callback instance, then printed on teardown. By printing on teardown, we
    do not interfere with the progress bar callback.
    """

    def __init__(self, **kwargs) -> None:
        warmup_epochs = kwargs.pop("warmup_epochs")  # a local parameter we added
        super().__init__(**kwargs)
        self.early_stopping_reason = None
        self.warmup_epochs = warmup_epochs
        self._current_step = 0  # internal counter of steps to check

    def _evaluate_stopping_criteria(self, current: torch.Tensor) -> tuple[bool, str]:
        self._current_step += 1
        if self._current_step <= self.warmup_epochs:
            # If we're in the warm-up phase, don't apply early stopping
            return False, "Early Stop not applied during first " + str(
                self.warmup_epochs
            ) + " steps"
        else:
            should_stop, reason = super()._evaluate_stopping_criteria(current)
            if should_stop:
                self.early_stopping_reason = reason
            return should_stop, reason

    def teardown(
        self,
        _trainer: pl.Trainer,
        _pl_module: pl.LightningModule,
        stage: str | None = None,
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


class ScibCallback(Callback):
    """A callback to initialize the Scib-Metrics autotune module."""

    def __init__(
        self,
    ):
        super().__init__()
        self.pl_module = None

    def _get_report_dict(self, trainer: pl.Trainer, pl_module: pl.LightningModule):
        """Exposing the pl_module to the scib-metrics autotune"""
        self.pl_module = pl_module

    def on_train_end(self, trainer: pl.Trainer, pl_module: pl.LightningModule):
        self._get_report_dict(trainer, pl_module)

    def on_validation_end(self, trainer: pl.Trainer, pl_module: pl.LightningModule):
        self._get_report_dict(trainer, pl_module)
