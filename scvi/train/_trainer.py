import sys
import warnings
from typing import Optional, Union

import numpy as np
import pytorch_lightning as pl
from pytorch_lightning.callbacks.early_stopping import EarlyStopping
from pytorch_lightning.loggers import LightningLoggerBase

from scvi import settings
from scvi._compat import Literal

from ._logger import SimpleLogger
from ._progress import ProgressBar
from ._trainingplans import PyroTrainingPlan


class Trainer(pl.Trainer):
    """
    Lightweight wrapper of Pytorch Lightning Trainer.

    Appropriate defaults are set for scvi-tools models, as well as callbacks like
    EarlyStopping, with parameters accessible here.

    Parameters
    ----------
    gpus
        Number of gpus to train on (int) or which GPUs to train on (list or str) applied per node
    benchmark
        If true enables cudnn.benchmark, which improves speed when inputs are fixed size
    flush_logs_every_n_steps
        How often to flush logs to disk. By default, flushes after training complete.
    check_val_every_n_epoch
        Check val every n train epochs. By default, val is not checked, unless `early_stopping` is `True`.
    max_epochs
        Stop training once this number of epochs is reached.
    default_root_dir
        Default path for logs and weights when no logger/ckpt_callback passed.
        Defaults to `scvi.settings.logging_dir`. Can be remote file paths such as
        s3://mybucket/path or ‘hdfs://path/’
    checkpoint_callback
        If `True`, enable checkpointing. It will configure a default ModelCheckpoint
        callback if there is no user-defined ModelCheckpoint in `callbacks`.
    num_sanity_val_steps
        Sanity check runs n validation batches before starting the training routine.
        Set it to -1 to run all batches in all validation dataloaders.
    weights_summary
        Prints a summary of the weights when training begins.
    early_stopping
        Whether to perform early stopping with respect to the validation set. This
        automatically adds a :class:`~pytorch_lightning.callbacks.early_stopping.EarlyStopping`
        instance. A custom instance can be passed by using the callbacks argument and
        setting this to `False`.
    early_stopping_monitor
        Metric logged during validation set epoch. The available metrics will depend on
        the training plan class used. We list the most common options here in the typing.
    early_stopping_min_delta
        Minimum change in the monitored quantity to qualify as an improvement,
        i.e. an absolute change of less than min_delta, will count as no improvement.
    early_stopping_patience
        Number of validation epochs with no improvement after which training will be stopped.
    early_stopping_mode
            In 'min' mode, training will stop when the quantity monitored has stopped decreasing
            and in 'max' mode it will stop when the quantity monitored has stopped increasing.
    progress_bar_refresh_rate
        How often to refresh progress bar (in steps). Value 0 disables progress bar.
    simple_progress_bar
        Use custom scvi-tools simple progress bar (per epoch rather than per batch)
    logger
        A valid pytorch lightning logger. Defaults to a simple dictionary logger.
        If `True`, defaults to the default pytorch lightning logger.
    **kwargs
        Other keyword args for :class:`~pytorch_lightning.trainer.Trainer`
    """

    def __init__(
        self,
        gpus: Union[int, str] = 1,
        benchmark: bool = True,
        flush_logs_every_n_steps=np.inf,
        check_val_every_n_epoch: Optional[int] = None,
        max_epochs: int = 400,
        default_root_dir: Optional[str] = None,
        checkpoint_callback: bool = False,
        num_sanity_val_steps: int = 0,
        weights_summary: Optional[Literal["top", "full"]] = None,
        early_stopping: bool = False,
        early_stopping_monitor: Literal[
            "elbo_validation", "reconstruction_loss_validation", "kl_local_validation"
        ] = "elbo_validation",
        early_stopping_min_delta: float = 0.00,
        early_stopping_patience: int = 45,
        early_stopping_mode: Literal["min", "max"] = "min",
        progress_bar_refresh_rate: int = 1,
        simple_progress_bar: bool = True,
        logger: Union[Optional[LightningLoggerBase], bool] = None,
        **kwargs
    ):
        if default_root_dir is None:
            default_root_dir = settings.logging_dir

        kwargs["callbacks"] = (
            [] if "callbacks" not in kwargs.keys() else kwargs["callbacks"]
        )
        if early_stopping:
            early_stopping_callback = EarlyStopping(
                monitor=early_stopping_monitor,
                min_delta=early_stopping_min_delta,
                patience=early_stopping_patience,
                mode=early_stopping_mode,
            )
            kwargs["callbacks"] += [early_stopping_callback]
            check_val_every_n_epoch = 1
        else:
            check_val_every_n_epoch = (
                check_val_every_n_epoch
                if check_val_every_n_epoch is not None
                # needs to be an integer, np.inf does not work
                else sys.maxsize
            )

        if simple_progress_bar:
            bar = ProgressBar(refresh_rate=progress_bar_refresh_rate)
            kwargs["callbacks"] += [bar]

        if logger is None:
            logger = SimpleLogger()

        super().__init__(
            gpus=gpus,
            benchmark=benchmark,
            flush_logs_every_n_steps=flush_logs_every_n_steps,
            check_val_every_n_epoch=check_val_every_n_epoch,
            max_epochs=max_epochs,
            default_root_dir=default_root_dir,
            checkpoint_callback=checkpoint_callback,
            num_sanity_val_steps=num_sanity_val_steps,
            weights_summary=weights_summary,
            logger=logger,
            progress_bar_refresh_rate=progress_bar_refresh_rate,
            **kwargs,
        )

    def fit(self, *args, **kwargs):

        with warnings.catch_warnings():
            warnings.filterwarnings(
                action="ignore", category=UserWarning, message="The dataloader,"
            )
            warnings.filterwarnings(
                action="ignore",
                category=UserWarning,
                message="you defined a validation_step but have no val_dataloader",
            )
            warnings.filterwarnings(
                action="ignore",
                category=UserWarning,
                message="One of given dataloaders is None and it will be skipped",
            )
            if isinstance(args[0], PyroTrainingPlan):
                warnings.filterwarnings(
                    action="ignore",
                    category=UserWarning,
                    message="`LightningModule.configure_optimizers` returned `None`",
                )
            super().fit(*args, **kwargs)
