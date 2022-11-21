import warnings
from typing import Optional, Union

import numpy as np

from scvi._decorators import classproperty
from scvi.dataloaders import DataSplitter
from scvi.train import TrainingPlan, TrainRunner


class BaseTrainingMixin:
    """Base training mixin class."""

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
        *,
        data_splitter_kwargs: Optional[dict] = None,
        training_plan_kwargs: Optional[dict] = None,
        train_runner_kwargs: Optional[dict] = None,
        catch_lightning_warnings: bool = False,
    ) -> None:
        """Train the model."""
        data_splitter_kwargs = data_splitter_kwargs or {}
        training_plan_kwargs = training_plan_kwargs or {}
        train_runner_kwargs = train_runner_kwargs or {}

        data_splitter = self._data_splitter_cls(
            self.adata_manager, **data_splitter_kwargs
        )
        training_plan = self._training_plan_cls(self.module, **training_plan_kwargs)
        train_runner = self._train_runner_cls(
            self,
            training_plan=training_plan,
            data_splitter=data_splitter,
            **train_runner_kwargs,
        )
        if catch_lightning_warnings:
            with warnings.catch_warnings():
                warnings.filterwarnings(
                    "ignore", category=UserWarning, module=r"pytorch_lightning.*"
                )
                train_runner()
        else:
            train_runner()


class UnsupervisedTrainingMixin(BaseTrainingMixin):
    """General purpose unsupervised train method."""

    def train(
        self,
        max_epochs: Optional[int] = None,
        use_gpu: Optional[Union[str, int, bool]] = None,
        train_size: float = 0.9,
        validation_size: Optional[float] = None,
        batch_size: int = 128,
        early_stopping: bool = False,
        plan_kwargs: Optional[dict] = None,
        **trainer_kwargs,
    ) -> None:
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
        plan_kwargs = plan_kwargs or {}
        trainer_kwargs = trainer_kwargs or {}

        if max_epochs is None:
            max_epochs = int(np.min([round((20000 / self.adata.n_obs) * 400), 400]))

        data_splitter_kwargs = {
            "train_size": train_size,
            "validation_size": validation_size,
            "batch_size": batch_size,
            "use_gpu": use_gpu,
        }
        trainer_kwargs["early_stopping"] = trainer_kwargs.get(
            "early_stopping", early_stopping
        )
        train_runner_kwargs = {
            "max_epochs": max_epochs,
            "use_gpu": use_gpu,
        }
        train_runner_kwargs.update(trainer_kwargs)
        super().train(
            data_splitter_kwargs=data_splitter_kwargs,
            training_plan_kwargs=plan_kwargs,
            train_runner_kwargs=train_runner_kwargs,
        )
