from __future__ import annotations

import logging
import traceback

import mlx.core as mx
import pandas as pd

from scvi import settings
from scvi.dataloaders import DataSplitter
from scvi.model._utils import get_max_epochs_heuristic
from scvi.train import MlxTrainingPlan, TrainRunner
from scvi.utils import is_package_installed, mlflow_logger
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

        if settings.mlflow_set_tracking_uri != "":
            if is_package_installed("mlflow"):
                import mlflow

                mlflow.set_tracking_uri(settings.mlflow_set_tracking_uri)
                mlflow.set_experiment(settings.mlflow_set_experiment)

                with mlflow.start_run(
                    run_name=getattr(self, "run_name", None), log_system_metrics=True
                ):
                    try:
                        self._run_training_loop(max_epochs, data_splitter)

                        self.run_id = mlflow.active_run().info.run_id

                        mlflow_logger(
                            model=self,
                            trainer=None,
                            training_plan=self.training_plan,
                            data_splitter=data_splitter,
                            run_id=self.run_id,
                        )
                        mlflow.log_params({"max_epochs": max_epochs}, run_id=self.run_id)

                        mlflow.set_tag("status", "success")

                    except Exception as e:
                        error_msg = "".join(
                            traceback.format_exception(type(e), e, e.__traceback__)
                        )
                        with open("error_log.txt", "w") as f:
                            f.write(error_msg)

                        mlflow.log_artifact("error_log.txt", artifact_path="errors")
                        mlflow.set_tag("status", "failed")
                        mlflow.log_param("error_type", type(e).__name__)
                        mlflow.log_param("error_message", str(e))
                        print(f"Error logged to MLflow: {e}")
                        raise
            else:
                raise ModuleNotFoundError("Please install mlflow to use this functionality.")
        else:
            self._run_training_loop(max_epochs, data_splitter)

        return self

    def _run_training_loop(self, max_epochs: int, data_splitter: DataSplitter):
        """Run the core MLX training loop.

        Parameters
        ----------
        max_epochs
            Number of training epochs.
        data_splitter
            Initialized DataSplitter with train/val dataloaders.
        """
        # Get data loaders
        train_dl = data_splitter.train_dataloader()
        val_dl = data_splitter.val_dataloader()

        # Training loop
        self.training_plan.current_epoch = 0
        self.training_plan.current_step = 0

        history = {
            "train_loss": [],
            "elbo_train": [],
            "reconstruction_loss_train": [],
            "kl_local_train": [],
            "kl_global_train": [],
        }

        for epoch in range(max_epochs):
            self.training_plan.current_epoch = epoch
            epoch_loss = 0.0
            epoch_rec_loss = 0.0
            epoch_kl_local = 0.0
            epoch_kl_global = 0.0
            epoch_n_obs = 0
            n_batches = 0

            # Training phase
            for batch in train_dl:
                self.training_plan.current_step += 1
                try:
                    output = self.training_plan.train_step(batch)
                    if mx.isnan(output["loss"]).any():
                        logger.warning("Skipping NaN batch")
                        continue
                    epoch_loss += output["loss"]
                    epoch_rec_loss += output.get("reconstruction_loss", 0.0)
                    epoch_kl_local += output.get("kl_local", 0.0)
                    epoch_kl_global += output.get("kl_global", 0.0)
                    epoch_n_obs += output.get("n_obs", 0)
                    n_batches += 1
                except Exception as e:
                    logger.error(f"Error processing training batch: {str(e)}")
                    raise

            avg_loss = epoch_loss / max(n_batches, 1)
            avg_loss_val = float(avg_loss)
            n_obs = max(epoch_n_obs, 1)
            history["train_loss"].append(avg_loss_val)
            history["elbo_train"].append(avg_loss_val)
            history["reconstruction_loss_train"].append(epoch_rec_loss / n_obs)
            history["kl_local_train"].append(epoch_kl_local / n_obs)
            history["kl_global_train"].append(epoch_kl_global / max(n_batches, 1))
            logger.info(f"Epoch {epoch + 1}/{max_epochs}, Loss: {avg_loss:.4f}")

            # Validation phase
            if val_dl is not None:
                val_loss = 0.0
                val_rec_loss = 0.0
                val_kl_local = 0.0
                val_kl_global = 0.0
                val_n_obs = 0
                n_val_batches = 0

                for batch in val_dl:
                    try:
                        output = self.training_plan.validate_step(batch)
                        val_loss += output["validation_loss"]
                        val_rec_loss += output.get("reconstruction_loss", 0.0)
                        val_kl_local += output.get("kl_local", 0.0)
                        val_kl_global += output.get("kl_global", 0.0)
                        val_n_obs += output.get("n_obs", 0)
                        n_val_batches += 1
                    except Exception as e:
                        logger.error(f"Error processing validation batch: {str(e)}")
                        raise

                avg_val_loss = val_loss / max(n_val_batches, 1)
                avg_val_loss_val = float(avg_val_loss)
                n_val = max(val_n_obs, 1)
                history.setdefault("validation_loss", []).append(avg_val_loss_val)
                history.setdefault("elbo_validation", []).append(avg_val_loss_val)
                history.setdefault("reconstruction_loss_validation", []).append(
                    val_rec_loss / n_val
                )
                history.setdefault("kl_local_validation", []).append(val_kl_local / n_val)
                history.setdefault("kl_global_validation", []).append(
                    val_kl_global / max(n_val_batches, 1)
                )
                logger.info(f"Validation loss: {avg_val_loss:.4f}")

        self.history_ = {k: pd.DataFrame({k: v}) for k, v in history.items()}

        # Store indices from data splitter
        self.train_indices = getattr(data_splitter, "train_idx", None)
        self.test_indices = getattr(data_splitter, "test_idx", None)
        self.validation_indices = getattr(data_splitter, "val_idx", None)

        # Set state after training completes
        self.is_trained_ = True
        if hasattr(self.module, "eval"):
            self.module.eval()
        else:
            self.training_plan._set_module_training(False)

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
