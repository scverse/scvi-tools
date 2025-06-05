import sys
import warnings
from typing import Literal

import lightning.pytorch as pl
from lightning.pytorch.accelerators import Accelerator
from lightning.pytorch.callbacks import LearningRateMonitor
from lightning.pytorch.loggers import Logger

from scvi import settings
from scvi.model._utils import use_distributed_sampler

from ._callbacks import (
    LoudEarlyStopping,
    SaveCheckpoint,
)
from ._logger import SimpleLogger
from ._progress import ProgressBar
from ._trainingplans import PyroTrainingPlan


class Trainer(pl.Trainer):
    """Lightweight wrapper of Pytorch Lightning Trainer.

    Appropriate defaults are set for scvi-tools models, as well as callbacks like
    EarlyStopping, with parameters accessible here.

    Parameters
    ----------
    accelerator
        Supports passing different accelerator types ("cpu", "gpu", "tpu", "ipu", "hpu", "mps,
        "auto") as well as custom accelerator instances.
    devices
        The devices to use. Can be set to a positive number (int or str), a sequence of device
        indices (list or str), the value ``-1`` to indicate all available devices should be used,
        or ``"auto"`` for automatic selection based on the chosen accelerator. Default: ``"auto"``.
    benchmark
        If true enables cudnn.benchmark, which improves speed when inputs are fixed size
    check_val_every_n_epoch
        Check val every n train epochs. By default, val is not checked, unless `early_stopping` is
        `True`.
    max_epochs
        Stop training once this number of epochs is reached.
    default_root_dir
        Default path for logs and weights when no logger/ckpt_callback passed.
        Defaults to `scvi.settings.logging_dir`. Can be remote file paths such as
        s3://mybucket/path or ‘hdfs://path/’
    enable_checkpointing
        If ``True``, enables checkpointing with a default :class:`~scvi.train.SaveCheckpoint`
        callback if there is no user-defined :class:`~scvi.train.SaveCheckpoint` in ``callbacks``.
    checkpointing_monitor
        If ``enable_checkpointing`` is ``True``, specifies the metric to monitor for checkpointing.
    num_sanity_val_steps
        Sanity check runs n validation batches before starting the training routine.
        Set it to -1 to run all batches in all validation dataloaders.
    enable_model_summary
        Whether to enable or disable the model summarization.
    early_stopping
        Whether to perform early stopping with respect to the validation set. This
        automatically adds a :class:`~lightning.pytorch.callbacks.early_stopping.EarlyStopping`
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
    early_stopping_warmup_epochs
        Wait for a certain number of warm-up epochs before the early stopping starts monitoring
    early_stopping_mode
        In 'min' mode, training will stop when the quantity monitored has stopped decreasing
        and in 'max' mode it will stop when the quantity monitored has stopped increasing.
    enable_progress_bar
        Whether to enable or disable the progress bar.
    progress_bar_refresh_rate
        How often to refresh progress bar (in steps). Value 0 disables progress bar.
    simple_progress_bar
        Use custom scvi-tools simple progress bar (per epoch rather than per batch).
        When `False`, uses default PyTorch Lightning progress bar, unless `enable_progress_bar`
        is `False`.
    logger
        A valid pytorch lightning logger. Defaults to a simple dictionary logger.
        If `True`, defaults to the default pytorch lightning logger.
    log_every_n_steps
        How often to log within steps. This does not affect epoch-level logging.
    **kwargs
        Other keyword args for :class:`~pytorch_lightning.trainer.Trainer`
    """

    def __init__(
        self,
        accelerator: str | Accelerator | None = None,
        devices: list[int] | str | int | None = None,
        benchmark: bool = True,
        check_val_every_n_epoch: int | None = None,
        max_epochs: int = 400,
        default_root_dir: str | None = None,
        enable_checkpointing: bool = False,
        checkpointing_monitor: str = "validation_loss",
        num_sanity_val_steps: int = 0,
        enable_model_summary: bool = False,
        early_stopping: bool = False,
        early_stopping_monitor: Literal[
            "elbo_validation", "reconstruction_loss_validation", "kl_local_validation"
        ] = "elbo_validation",
        early_stopping_min_delta: float = 0.00,
        early_stopping_patience: int = 45,
        early_stopping_warmup_epochs: int = 0,
        early_stopping_mode: Literal["min", "max"] = "min",
        enable_progress_bar: bool = True,
        progress_bar_refresh_rate: int = 1,
        simple_progress_bar: bool = True,
        logger: Logger | None | bool = None,
        log_every_n_steps: int = 10,
        learning_rate_monitor: bool = False,
        **kwargs,
    ):
        if default_root_dir is None:
            default_root_dir = settings.logging_dir

        check_val_every_n_epoch = check_val_every_n_epoch or sys.maxsize
        callbacks = kwargs.pop("callbacks", [])

        if use_distributed_sampler(kwargs.get("strategy", None)):
            warnings.warn(
                "early_stopping was automaticaly disabled due to the use of DDP",
                UserWarning,
                stacklevel=settings.warnings_stacklevel,
            )
            early_stopping = False

        if early_stopping:
            early_stopping_callback = LoudEarlyStopping(
                monitor=early_stopping_monitor,
                min_delta=early_stopping_min_delta,
                patience=early_stopping_patience,
                mode=early_stopping_mode,
                warmup_epochs=early_stopping_warmup_epochs,
            )
            callbacks.append(early_stopping_callback)
            check_val_every_n_epoch = 1

        if enable_checkpointing and not any(isinstance(c, SaveCheckpoint) for c in callbacks):
            callbacks.append(SaveCheckpoint(monitor=checkpointing_monitor))
            check_val_every_n_epoch = 1
        elif any(isinstance(c, SaveCheckpoint) for c in callbacks):
            # check if user provided already provided the callback
            enable_checkpointing = True
            check_val_every_n_epoch = 1

        if learning_rate_monitor and not any(
            isinstance(c, LearningRateMonitor) for c in callbacks
        ):
            callbacks.append(LearningRateMonitor())
            check_val_every_n_epoch = 1

        if simple_progress_bar and enable_progress_bar:
            callbacks.append(ProgressBar(refresh_rate=progress_bar_refresh_rate))

        if logger is None:
            logger = SimpleLogger()

        super().__init__(
            accelerator=accelerator,
            devices=devices,
            benchmark=benchmark,
            check_val_every_n_epoch=check_val_every_n_epoch,
            max_epochs=max_epochs,
            default_root_dir=default_root_dir,
            enable_checkpointing=enable_checkpointing,
            num_sanity_val_steps=num_sanity_val_steps,
            enable_model_summary=enable_model_summary,
            logger=logger,
            log_every_n_steps=log_every_n_steps,
            enable_progress_bar=enable_progress_bar,
            callbacks=callbacks,
            **kwargs,
        )

    def fit(self, *args, **kwargs):
        """Fit the model."""
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
            # bug in pytorch lightning, assumes SequentialSampler
            # https://github.com/PyTorchLightning/pytorch-lightning/blob/
            # 48cb38ac5dd0159c8f7c5189c888dfd04a2ed34b/pytorch_lightning/
            # trainer/data_loading.py#L311-L314
            warnings.filterwarnings(
                action="ignore",
                category=UserWarning,
                message="Your `val_dataloader` has `shuffle=True`",
            )
            if isinstance(args[0], PyroTrainingPlan):
                warnings.filterwarnings(
                    action="ignore",
                    category=UserWarning,
                    message="`LightningModule.configure_optimizers` returned `None`",
                )
            try:
                super().fit(*args, **kwargs)
            except NameError:
                import gc

                gc.collect()
                import torch

                if torch.cuda.is_available():
                    torch.cuda.empty_cache()
