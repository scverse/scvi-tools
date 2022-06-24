import warnings
from asyncio.log import logger
from typing import Optional

import jax
import numpy as np

from scvi.dataloaders import DataSplitter
from scvi.train import JaxModuleInit, JaxTrainingPlan, TrainRunner


class JaxTrainingMixin:
    def train(
        self,
        max_epochs: Optional[int] = None,
        use_gpu: Optional[bool] = None,
        train_size: float = 0.9,
        validation_size: Optional[float] = None,
        batch_size: int = 128,
        lr: float = 1e-3,
        **trainer_kwargs,
    ):
        if max_epochs is None:
            n_cells = self.adata.n_obs
            max_epochs = np.min([round((20000 / n_cells) * 400), 400])

        if use_gpu is None or use_gpu is True:
            try:
                self.module.to(jax.devices("gpu")[0])
                logger.debug(
                    "Jax module moved to GPU. "
                    "Note: Pytorch lightning will show GPU is not being used for the Trainer."
                )
            except RuntimeError:
                logger.debug("No GPU available to Jax.")
        else:
            cpu_device = jax.devices("cpu")[0]
            if self.module.device is not cpu_device:
                self.module.to(cpu_device)
                logger.debug("Jax module moved back to CPU. ")

        data_splitter = DataSplitter(
            self.adata_manager,
            train_size=train_size,
            validation_size=validation_size,
            batch_size=batch_size,
            # for pinning memory only
            use_gpu=False,
            iter_ndarray=True,
        )

        self.training_plan = JaxTrainingPlan(
            self.module, optim_kwargs=dict(learning_rate=lr)
        )
        if "callbacks" not in trainer_kwargs.keys():
            trainer_kwargs["callbacks"] = []
        trainer_kwargs["callbacks"].append(JaxModuleInit())

        # Ignore Pytorch Lightning warnings for Jax workarounds.
        with warnings.catch_warnings():
            warnings.filterwarnings(
                "ignore", category=UserWarning, module=r"pytorch_lightning.*"
            )
            runner = TrainRunner(
                self,
                training_plan=self.training_plan,
                data_splitter=data_splitter,
                max_epochs=max_epochs,
                use_gpu=use_gpu,
                **trainer_kwargs,
            )
            runner()

        self.is_trained_ = True
        self.module.eval()
