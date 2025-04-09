from __future__ import annotations

import logging
import warnings

from scvi.dataloaders import DataSplitter
from scvi.model._utils import get_max_epochs_heuristic, parse_device_args
from scvi.train import JaxModuleInit, JaxTrainingPlan, TrainRunner
from scvi.utils._docstrings import devices_dsp

logger = logging.getLogger(__name__)


class JaxTrainingMixin:
    """General purpose train method for Jax-backed modules."""

    _data_splitter_cls = DataSplitter
    _training_plan_cls = JaxTrainingPlan
    _train_runner_cls = TrainRunner

    @devices_dsp.dedent
    def train(
        self,
        max_epochs: int | None = None,
        accelerator: str = "auto",
        devices: int | list[int] | str = "auto",
        train_size: float | None = None,
        validation_size: float | None = None,
        shuffle_set_split: bool = True,
        batch_size: int = 128,
        datasplitter_kwargs: dict | None = None,
        plan_kwargs: dict | None = None,
        **trainer_kwargs,
    ):
        """Train the model.

        Parameters
        ----------
        max_epochs
            Number of passes through the dataset. If `None`, defaults to
            `np.min([round((20000 / n_cells) * 400), 400])`
        %(param_accelerator)s
        %(param_devices)s
        train_size
            Size of training set in the range [0.0, 1.0].
        validation_size
            Size of the test set. If `None`, defaults to 1 - `train_size`. If
            `train_size + validation_size < 1`, the remaining cells belong to a test set.
        shuffle_set_split
            Whether to shuffle indices before splitting. If `False`, the val, train, and test set
            are split in the sequential order of the data according to `validation_size` and
            `train_size` percentages.
        batch_size
            Minibatch size to use during training.
        lr
            Learning rate to use during training.
        datasplitter_kwargs
            Additional keyword arguments passed into :class:`~scvi.dataloaders.DataSplitter`.
        plan_kwargs
            Keyword args for :class:`~scvi.train.JaxTrainingPlan`. Keyword arguments passed to
            `train()` will overwrite values present in `plan_kwargs`, when appropriate.
        **trainer_kwargs
            Other keyword args for :class:`~scvi.train.Trainer`.
        """
        if max_epochs is None:
            max_epochs = get_max_epochs_heuristic(self.adata.n_obs)

        _, _, device = parse_device_args(
            accelerator,
            devices,
            return_device="jax",
            validate_single_device=True,
        )
        try:
            self.module.to(device)
            logger.info(
                f"Jax module moved to {device}."
                "Note: Pytorch lightning will show GPU is not being used for the Trainer."
            )
        except RuntimeError:
            logger.debug("No GPU available to Jax.")

        datasplitter_kwargs = datasplitter_kwargs or {}

        data_splitter = self._data_splitter_cls(
            self.adata_manager,
            train_size=train_size,
            validation_size=validation_size,
            shuffle_set_split=shuffle_set_split,
            batch_size=batch_size,
            iter_ndarray=True,
            **datasplitter_kwargs,
        )
        plan_kwargs = plan_kwargs if isinstance(plan_kwargs, dict) else {}

        self.training_plan = self._training_plan_cls(self.module, **plan_kwargs)
        if "callbacks" not in trainer_kwargs.keys():
            trainer_kwargs["callbacks"] = []
        trainer_kwargs["callbacks"].append(JaxModuleInit())

        # Ignore Pytorch Lightning warnings for Jax workarounds.
        with warnings.catch_warnings():
            warnings.filterwarnings("ignore", category=UserWarning, module=r"pytorch_lightning.*")
            runner = self._train_runner_cls(
                self,
                training_plan=self.training_plan,
                data_splitter=data_splitter,
                max_epochs=max_epochs,
                accelerator="cpu",
                devices="auto",
                **trainer_kwargs,
            )
            runner()

        self.is_trained_ = True
        self.module.eval()
