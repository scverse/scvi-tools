from __future__ import annotations

import logging
import warnings
from typing import TYPE_CHECKING

from scvi import settings
from scvi.dataloaders import DataSplitter
from scvi.model._utils import get_max_epochs_heuristic, parse_device_args
from scvi.train import JaxModuleInit, JaxTrainingPlan, TrainRunner
from scvi.train._config import merge_kwargs
from scvi.utils._docstrings import devices_dsp

if TYPE_CHECKING:
    from scvi.train._config import KwargsLike

logger = logging.getLogger(__name__)


class JaxTrainingMixin:
    """General purpose train method for Jax-backed modules."""

    _data_splitter_cls = DataSplitter
    _training_plan_cls = JaxTrainingPlan
    _train_runner_cls = TrainRunner

    @staticmethod
    def _resolve_jax_n_devices(accelerator: str, devices) -> int:
        """Resolve the number of JAX devices to use based on user arguments."""
        import jax

        if accelerator == "cpu":
            return 1
        available = jax.local_device_count()
        if devices == "auto" or devices == -1:
            return available
        if isinstance(devices, int):
            return min(devices, available)
        if isinstance(devices, (list, tuple)):
            return min(len(devices), available)
        if isinstance(devices, str) and devices.isdigit():
            return min(int(devices), available)
        return 1

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
        plan_config: KwargsLike | None = None,
        plan_kwargs: KwargsLike | None = None,
        trainer_config: KwargsLike | None = None,
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
        plan_config
            Configuration object or mapping used to build :class:`~scvi.train.JaxTrainingPlan`.
            Values in ``plan_kwargs`` and explicit arguments take precedence.
        trainer_config
            Configuration object or mapping used to build :class:`~scvi.train.Trainer`. Values in
            ``trainer_kwargs`` and explicit arguments take precedence.
        **trainer_kwargs
            Other keyword args for :class:`~scvi.train.Trainer`.
        """
        if max_epochs is None:
            max_epochs = get_max_epochs_heuristic(self.adata.n_obs)

        n_devices = self._resolve_jax_n_devices(accelerator, devices)

        _, _, device = parse_device_args(
            accelerator,
            devices,
            return_device="jax",
            validate_single_device=False,
        )
        try:
            self.module.to(device)
            logger.info(
                f"Jax module moved to {device}."
                "Note: Pytorch lightning will show GPU is not being used for the Trainer."
            )
        except RuntimeError:
            logger.debug("No GPU available to Jax.")

        if n_devices > 1:
            logger.info(f"JAX multi-GPU training with {n_devices} devices.")

        datasplitter_kwargs = datasplitter_kwargs or {}

        # For multi-GPU, ensure batch size is divisible by n_devices
        effective_batch_size = batch_size or settings.batch_size
        if n_devices > 1 and effective_batch_size % n_devices != 0:
            effective_batch_size = (effective_batch_size // n_devices + 1) * n_devices
            logger.info(
                f"Adjusted batch size to {effective_batch_size} for even sharding "
                f"across {n_devices} devices."
            )

        # For multi-GPU, drop incomplete training batches to avoid sharding issues
        if n_devices > 1:
            datasplitter_kwargs.setdefault("drop_last", True)

        data_splitter = self._data_splitter_cls(
            self.adata_manager,
            train_size=train_size,
            validation_size=validation_size,
            shuffle_set_split=shuffle_set_split,
            batch_size=effective_batch_size,
            iter_ndarray=True,
            **datasplitter_kwargs,
        )
        plan_kwargs = merge_kwargs(plan_config, plan_kwargs, name="plan")

        self.training_plan = self._training_plan_cls(self.module, **plan_kwargs)
        self.training_plan.n_devices = n_devices
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
                trainer_config=trainer_config,
                **trainer_kwargs,
            )
            runner()

        # After training, unreplicate state if multi-GPU was used
        if n_devices > 1 and self.module.train_state is not None:
            from flax.jax_utils import unreplicate

            self.module.train_state = unreplicate(self.module.train_state)
            self.training_plan.n_devices = 1

        self.is_trained_ = True
        self.module.eval()
