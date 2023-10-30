import logging
import warnings
from typing import Union

import lightning.pytorch as pl
import numpy as np
import pandas as pd

from scvi import settings
from scvi.dataloaders import DataSplitter, SemiSupervisedDataSplitter
from scvi.model._utils import parse_device_args
from scvi.model.base import BaseModelClass
from scvi.train import Trainer

logger = logging.getLogger(__name__)


class TrainRunner:
    """TrainRunner calls Trainer.fit() and handles pre and post training procedures.

    Parameters
    ----------
    model
        model to train
    training_plan
        initialized TrainingPlan
    data_splitter
        initialized :class:`~scvi.dataloaders.SemiSupervisedDataSplitter` or
        :class:`~scvi.dataloaders.DataSplitter`
    max_epochs
        max_epochs to train for
    accelerator
        Supports passing different accelerator types ("cpu", "gpu", "tpu", "ipu", "hpu",
        "mps, "auto") as well as custom accelerator instances.
    devices
        The devices to use. Can be set to a positive number (int or str), a sequence of
        device indices (list or str), the value -1 to indicate all available devices should
        be used, or "auto" for automatic selection based on the chosen accelerator.
    trainer_kwargs
        Extra kwargs for :class:`~scvi.train.Trainer`

    Examples
    --------
    >>> # Following code should be within a subclass of BaseModelClass
    >>> data_splitter = DataSplitter(self.adata)
    >>> training_plan = TrainingPlan(self.module, len(data_splitter.train_idx))
    >>> runner = TrainRunner(
    >>>     self,
    >>>     training_plan=trianing_plan,
    >>>     data_splitter=data_splitter,
    >>>     max_epochs=max_epochs)
    >>> runner()
    """

    _trainer_cls = Trainer

    def __init__(
        self,
        model: BaseModelClass,
        training_plan: pl.LightningModule,
        data_splitter: Union[SemiSupervisedDataSplitter, DataSplitter],
        max_epochs: int,
        accelerator: str = "auto",
        devices: Union[int, list[int], str] = "auto",
        **trainer_kwargs,
    ):
        self.training_plan = training_plan
        self.data_splitter = data_splitter
        self.model = model
        accelerator, lightning_devices, device = parse_device_args(
            accelerator=accelerator,
            devices=devices,
            return_device="torch",
        )
        self.accelerator = accelerator
        self.lightning_devices = lightning_devices
        self.device = device

        if getattr(self.training_plan, "reduce_lr_on_plateau", False):
            trainer_kwargs["learning_rate_monitor"] = True

        self.trainer = self._trainer_cls(
            max_epochs=max_epochs,
            accelerator=accelerator,
            devices=lightning_devices,
            **trainer_kwargs,
        )
        # currently set for MetricsCallback
        self.trainer._model = model  # TODO: Find a better way to do this

    def __call__(self):
        """Run training."""
        if hasattr(self.data_splitter, "n_train"):
            self.training_plan.n_obs_training = self.data_splitter.n_train
        if hasattr(self.data_splitter, "n_val"):
            self.training_plan.n_obs_validation = self.data_splitter.n_val

        self.trainer.fit(self.training_plan, self.data_splitter)
        self._update_history()

        # data splitter only gets these attrs after fit
        self.model.train_indices = self.data_splitter.train_idx
        self.model.test_indices = self.data_splitter.test_idx
        self.model.validation_indices = self.data_splitter.val_idx

        self.model.module.eval()
        self.model.is_trained_ = True
        self.model.to_device(self.device)
        self.model.trainer = self.trainer

    def _update_history(self):
        # model is being further trained
        # this was set to true during first training session
        if self.model.is_trained_ is True:
            # if not using the default logger (e.g., tensorboard)
            if not isinstance(self.model.history_, dict):
                warnings.warn(
                    "Training history cannot be updated. Logger can be accessed from "
                    "`model.trainer.logger`",
                    UserWarning,
                    stacklevel=settings.warnings_stacklevel,
                )
                return
            else:
                new_history = self.trainer.logger.history
                for key, val in self.model.history_.items():
                    # e.g., no validation loss due to training params
                    if key not in new_history:
                        continue
                    prev_len = len(val)
                    new_len = len(new_history[key])
                    index = np.arange(prev_len, prev_len + new_len)
                    new_history[key].index = index
                    self.model.history_[key] = pd.concat(
                        [
                            val,
                            new_history[key],
                        ]
                    )
                    self.model.history_[key].index.name = val.index.name
        else:
            # set history_ attribute if it exists
            # other pytorch lightning loggers might not have history attr
            try:
                self.model.history_ = self.trainer.logger.history
            except AttributeError:
                self.history_ = None
