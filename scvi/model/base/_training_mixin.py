import warnings
from math import ceil
from typing import Any, Dict, Optional, Union

import numpy as np

from scvi.dataloaders import DataSplitter
from scvi.train import TrainingPlan, TrainRunner


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
            max_epochs = np.min([round((20000 / n_cells) * 400), 400])

        plan_kwargs = plan_kwargs if isinstance(plan_kwargs, dict) else dict()

        _check_warmup(plan_kwargs, max_epochs, n_cells, batch_size)

        data_splitter = DataSplitter(
            self.adata_manager,
            train_size=train_size,
            validation_size=validation_size,
            batch_size=batch_size,
            use_gpu=use_gpu,
        )
        training_plan = TrainingPlan(self.module, **plan_kwargs)

        es = "early_stopping"
        trainer_kwargs[es] = (
            early_stopping if es not in trainer_kwargs.keys() else trainer_kwargs[es]
        )
        runner = TrainRunner(
            self,
            training_plan=training_plan,
            data_splitter=data_splitter,
            max_epochs=max_epochs,
            use_gpu=use_gpu,
            **trainer_kwargs,
        )
        return runner()
