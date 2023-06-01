import logging
import warnings
from typing import List, Optional, Union

import jax

from scvi.dataloaders import DataSplitter
from scvi.model._utils import get_max_epochs_heuristic
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
        max_epochs: Optional[int] = None,
        use_gpu: Optional[Union[str, int, bool]] = None,
        accelerator: str = "auto",
        devices: Union[int, List[int], str] = "auto",
        train_size: float = 0.9,
        validation_size: Optional[float] = None,
        shuffle_set_split: bool = True,
        batch_size: int = 128,
        plan_kwargs: Optional[dict] = None,
        **trainer_kwargs,
    ):
        """Train the model.

        Parameters
        ----------
        max_epochs
            Number of passes through the dataset. If `None`, defaults to
            `np.min([round((20000 / n_cells) * 400), 400])`
        %(param_use_gpu)s
        %(param_accelerator)s
        %(param_devices)s
        train_size
            Size of training set in the range [0.0, 1.0].
        validation_size
            Size of the test set. If `None`, defaults to 1 - `train_size`. If
            `train_size + validation_size < 1`, the remaining cells belong to a test set.
        shuffle_set_split
            Whether to shuffle indices before splitting. If `False`, the val, train, and test set are split in the
            sequential order of the data according to `validation_size` and `train_size` percentages.
        batch_size
            Minibatch size to use during training.
        lr
            Learning rate to use during training.
        plan_kwargs
            Keyword args for :class:`~scvi.train.JaxTrainingPlan`. Keyword arguments passed to
            `train()` will overwrite values present in `plan_kwargs`, when appropriate.
        **trainer_kwargs
            Other keyword args for :class:`~scvi.train.Trainer`.
        """
        if max_epochs is None:
            max_epochs = get_max_epochs_heuristic(self.adata.n_obs)

        if use_gpu is None or use_gpu is True:
            try:
                self.module.to(jax.devices("gpu")[0])
                logger.info(
                    "Jax module moved to GPU. "
                    "Note: Pytorch lightning will show GPU is not being used for the Trainer."
                )
            except RuntimeError:
                logger.debug("No GPU available to Jax.")
        else:
            cpu_device = jax.devices("cpu")[0]
            self.module.to(cpu_device)
            logger.info("Jax module moved to CPU.")

        data_splitter = self._data_splitter_cls(
            self.adata_manager,
            train_size=train_size,
            validation_size=validation_size,
            shuffle_set_split=shuffle_set_split,
            batch_size=batch_size,
            iter_ndarray=True,
        )
        plan_kwargs = plan_kwargs if isinstance(plan_kwargs, dict) else {}

        self.training_plan = self._training_plan_cls(self.module, **plan_kwargs)
        if "callbacks" not in trainer_kwargs.keys():
            trainer_kwargs["callbacks"] = []
        trainer_kwargs["callbacks"].append(JaxModuleInit())

        # Ignore Pytorch Lightning warnings for Jax workarounds.
        with warnings.catch_warnings():
            warnings.filterwarnings(
                "ignore", category=UserWarning, module=r"pytorch_lightning.*"
            )
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
