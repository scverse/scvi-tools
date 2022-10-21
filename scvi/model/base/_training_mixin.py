import logging
import warnings
from math import ceil
from typing import Any, Dict, Optional, Union

import numpy as np

from scvi._decorators import classproperty
from scvi.autotune._types import Tunable
from scvi.dataloaders import DataSplitter, SemiSupervisedDataSplitter
from scvi.train import SemiSupervisedTrainingPlan, TrainingPlan, TrainRunner
from scvi.train._callbacks import SubSampleLabels

logger = logging.getLogger(__name__)


def _check_warmup(
    plan_kwargs: Dict[str, Any],
    max_epochs: int,
    n_cells: int,
    batch_size: int,
    train_size: float = 1.0,
) -> None:
    """
    Raises a warning if the max_kl_weight is not reached by the end of training.

    Parameters
    ----------
    plan_kwargs
        Keyword args for :class:`~scvi.train.TrainingPlan`.
    max_epochs
        Number of passes through the dataset.
    n_cells
        Number of cells in the whole datasets.
    batch_size
        Minibatch size to use during training.
    train_size
        Fraction of cells used for training.
    """
    _WARNING_MESSAGE = (
        "max_{mode}={max} is less than n_{mode}_kl_warmup={warm_up}. "
        "The max_kl_weight will not be reached during training."
    )

    n_steps_kl_warmup = plan_kwargs.get("n_steps_kl_warmup", None)
    n_epochs_kl_warmup = plan_kwargs.get("n_epochs_kl_warmup", None)

    # The only time n_steps_kl_warmup is used is when n_epochs_kl_warmup is explicitly
    # set to None. This also catches the case when both n_epochs_kl_warmup and
    # n_steps_kl_warmup are set to None and max_kl_weight will always be reached.
    if (
        "n_epochs_kl_warmup" in plan_kwargs
        and plan_kwargs["n_epochs_kl_warmup"] is None
    ):
        n_cell_train = ceil(train_size * n_cells)
        steps_per_epoch = n_cell_train // batch_size + (n_cell_train % batch_size >= 3)
        max_steps = max_epochs * steps_per_epoch
        if n_steps_kl_warmup and max_steps < n_steps_kl_warmup:
            warnings.warn(
                _WARNING_MESSAGE.format(
                    mode="steps", max=max_steps, warm_up=n_steps_kl_warmup
                )
            )
    elif n_epochs_kl_warmup:
        if max_epochs < n_epochs_kl_warmup:
            warnings.warn(
                _WARNING_MESSAGE.format(
                    mode="epochs", max=max_epochs, warm_up=n_epochs_kl_warmup
                )
            )
    else:
        if max_epochs < 400:
            warnings.warn(
                _WARNING_MESSAGE.format(mode="epochs", max=max_epochs, warm_up=400)
            )


class UnsupervisedTrainingMixin:
    """General purpose unsupervised train method."""

    @classproperty
    def _data_splitter_cls(cls) -> DataSplitter:
        return DataSplitter

    @classproperty
    def _training_plan_cls(cls) -> TrainingPlan:
        return TrainingPlan

    @classproperty
    def _train_runner_cls(cls) -> TrainRunner:
        return TrainRunner

    def train(
        self,
        max_epochs: Optional[int] = None,
        use_gpu: Optional[Union[str, int, bool]] = None,
        train_size: float = 0.9,
        validation_size: Optional[float] = None,
        batch_size: Tunable[int] = 128,
        early_stopping: bool = False,
        plan_kwargs: Optional[dict] = None,
        **trainer_kwargs,
    ):
        """
        Train the model.

        Parameters
        ----------
        max_epochs
            Number of passes through the dataset. If `None`, defaults to
            `np.min([round((20000 / n_cells) * 400), 400])`
        use_gpu
            Use default GPU if available (if None or True), or index of GPU to use (if int),
            or name of GPU (if str, e.g., `'cuda:0'`), or use CPU (if False).
        train_size
            Size of training set in the range [0.0, 1.0].
        validation_size
            Size of the test set. If `None`, defaults to 1 - `train_size`. If
            `train_size + validation_size < 1`, the remaining cells belong to a test set.
        batch_size
            Minibatch size to use during training.
        early_stopping
            Perform early stopping. Additional arguments can be passed in `**kwargs`.
            See :class:`~scvi.train.Trainer` for further options.
        plan_kwargs
            Keyword args for :class:`~scvi.train.TrainingPlan`. Keyword arguments passed to
            `train()` will overwrite values present in `plan_kwargs`, when appropriate.
        **trainer_kwargs
            Other keyword args for :class:`~scvi.train.Trainer`.
        """
        n_cells = self.adata.n_obs
        if max_epochs is None:
            max_epochs = int(np.min([round((20000 / n_cells) * 400), 400]))

        plan_kwargs = plan_kwargs if isinstance(plan_kwargs, dict) else dict()

        _check_warmup(plan_kwargs, max_epochs, n_cells, batch_size)

        data_splitter = self._data_splitter_cls(
            self.adata_manager,
            train_size=train_size,
            validation_size=validation_size,
            batch_size=batch_size,
            use_gpu=use_gpu,
        )
        training_plan = self._training_plan_cls(self.module, **plan_kwargs)

        es = "early_stopping"
        trainer_kwargs[es] = (
            early_stopping if es not in trainer_kwargs.keys() else trainer_kwargs[es]
        )
        runner = self._train_runner_cls(
            self,
            training_plan=training_plan,
            data_splitter=data_splitter,
            max_epochs=max_epochs,
            use_gpu=use_gpu,
            **trainer_kwargs,
        )
        return runner()


class SemiSupervisedTrainingMixin:
    """General purpose semi-supervised train method."""

    @classproperty
    def _data_splitter_cls(cls) -> DataSplitter:
        return SemiSupervisedDataSplitter

    @classproperty
    def _training_plan_cls(cls) -> TrainingPlan:
        return SemiSupervisedTrainingPlan

    @classproperty
    def _train_runner_cls(cls) -> TrainRunner:
        return TrainRunner

    def train(
        self,
        max_epochs: Optional[int] = None,
        n_samples_per_label: Optional[float] = None,
        check_val_every_n_epoch: Optional[int] = None,
        train_size: float = 0.9,
        validation_size: Optional[float] = None,
        batch_size: int = 128,
        use_gpu: Optional[Union[str, int, bool]] = None,
        plan_kwargs: Optional[dict] = None,
        **trainer_kwargs,
    ):
        """
        Train the model.

        Parameters
        ----------
        max_epochs
            Number of passes through the dataset for semisupervised training.
        n_samples_per_label
            Number of subsamples for each label class to sample per epoch. By default, there
            is no label subsampling.
        check_val_every_n_epoch
            Frequency with which metrics are computed on the data for validation set for both
            the unsupervised and semisupervised trainers. If you'd like a different frequency for
            the semisupervised trainer, set check_val_every_n_epoch in semisupervised_train_kwargs.
        train_size
            Size of training set in the range [0.0, 1.0].
        validation_size
            Size of the test set. If `None`, defaults to 1 - `train_size`. If
            `train_size + validation_size < 1`, the remaining cells belong to a test set.
        batch_size
            Minibatch size to use during training.
        use_gpu
            Use default GPU if available (if None or True), or index of GPU to use (if int),
            or name of GPU (if str, e.g., `'cuda:0'`), or use CPU (if False).
        plan_kwargs
            Keyword args for :class:`~scvi.train.SemiSupervisedTrainingPlan`. Keyword arguments passed to
            `train()` will overwrite values present in `plan_kwargs`, when appropriate.
        **trainer_kwargs
            Other keyword args for :class:`~scvi.train.Trainer`.
        """
        if max_epochs is None:
            n_cells = self.adata.n_obs
            max_epochs = int(np.min([round((20000 / n_cells) * 400), 400]))

            if self.was_pretrained:
                max_epochs = int(np.min([10, np.max([2, round(max_epochs / 3.0)])]))

        logger.info(f"Training for {max_epochs} epochs.")

        plan_kwargs = {} if plan_kwargs is None else plan_kwargs

        # if we have labeled cells, we want to subsample labels each epoch
        sampler_callback = (
            [SubSampleLabels()] if len(self._labeled_indices) != 0 else []
        )

        data_splitter = self._data_splitter_cls(
            adata_manager=self.adata_manager,
            train_size=train_size,
            validation_size=validation_size,
            n_samples_per_label=n_samples_per_label,
            batch_size=batch_size,
            use_gpu=use_gpu,
        )
        training_plan = self._training_plan_cls(self.module, **plan_kwargs)
        if "callbacks" in trainer_kwargs.keys():
            trainer_kwargs["callbacks"].concatenate(sampler_callback)
        else:
            trainer_kwargs["callbacks"] = sampler_callback

        runner = self._train_runner_cls(
            self,
            training_plan=training_plan,
            data_splitter=data_splitter,
            max_epochs=max_epochs,
            use_gpu=use_gpu,
            check_val_every_n_epoch=check_val_every_n_epoch,
            **trainer_kwargs,
        )
        return runner()
