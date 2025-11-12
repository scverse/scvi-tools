from __future__ import annotations

import logging

from scvi.dataloaders import DataSplitter
from scvi.model._utils import get_max_epochs_heuristic
from scvi.train import MlxTrainingPlan, TrainRunner
from scvi.utils._docstrings import devices_dsp

logger = logging.getLogger(__name__)


class MlxTrainingMixin:
    """Training mixin class for MLX models.

    This mixin class provides training functionality for models using MLX as a backend.
    It handles data splitting, training loops, and evaluation.
    """

    _data_splitter_cls = DataSplitter
    _training_plan_cls = MlxTrainingPlan
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
            Number of passes through the dataset. If None, defaults to
            `np.min([round((20000 / n_cells) * 400), 400])`
        accelerator
            Accelerator type. MLX automatically selects the best available device.
        devices
            Device selection. MLX automatically selects the best available device.
        train_size
            Training set size in the range [0.0, 1.0].
        validation_size
            Validation set size. If None, defaults to 1 - train_size.
        shuffle_set_split
            Whether to shuffle indices before splitting.
        batch_size
            Minibatch size to use during training.
        datasplitter_kwargs
            Additional keyword arguments for DataSplitter.
        plan_kwargs
            Keyword arguments for the training plan.
        **trainer_kwargs
            Additional keyword arguments for training.

        Returns
        -------
        self
            The trained model instance.
        """
        if max_epochs is None:
            max_epochs = get_max_epochs_heuristic(self.adata.n_obs)

        logger.info(f"Training model using MLX for {max_epochs} epochs")

        # Create data splitter
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

        # Setup data splitter
        data_splitter.setup()

        # Create the training plan
        plan_kwargs = plan_kwargs if isinstance(plan_kwargs, dict) else {}
        self.training_plan = self._training_plan_cls(self.module, **plan_kwargs)

        # Get data loaders
        train_dl = data_splitter.train_dataloader()
        val_dl = data_splitter.val_dataloader()

        # Training loop
        self.training_plan.current_epoch = 0
        self.training_plan.current_step = 0

        for epoch in range(max_epochs):
            self.training_plan.current_epoch = epoch
            epoch_loss = 0.0
            n_batches = 0

            # Training phase
            for batch in train_dl:
                self.training_plan.current_step += 1
                try:
                    output = self.training_plan.train_step(batch)
                    epoch_loss += output["loss"]
                    n_batches += 1
                except Exception as e:  # noqa: BLE001
                    logger.error(f"Error processing training batch: {str(e)}")
                    continue

            avg_loss = epoch_loss / max(n_batches, 1)  # Avoid division by zero
            logger.info(f"Epoch {epoch + 1}/{max_epochs}, Loss: {avg_loss:.4f}")

            # Validation phase
            if val_dl is not None:
                val_loss = 0.0
                n_val_batches = 0

                for batch in val_dl:
                    try:
                        output = self.training_plan.validate_step(batch)
                        val_loss += output["validation_loss"]
                        n_val_batches += 1
                    except Exception as e:  # noqa: BLE001
                        logger.error(f"Error processing validation batch: {str(e)}")
                        continue

                avg_val_loss = val_loss / max(n_val_batches, 1)  # Avoid division by zero
                logger.info(f"Validation loss: {avg_val_loss:.4f}")

        # Set state after training completes
        self.is_trained_ = True
        if hasattr(self.module, "eval"):
            self.module.eval()
        else:
            self.training_plan._set_module_training(False)

        return self

    def to_device(self, device):
        """Move the model to a specific device.

        MLX automatically handles device placement, so this is a no-op.

        Parameters
        ----------
        device
            Target device.
        """
        logger.info("MLX automatically handles device placement, ignoring to_device call")
        pass

    @property
    def device(self):
        """Get the current device.

        MLX automatically handles device placement.

        Returns
        -------
        str
            Device identifier.
        """
        return "mlx"
