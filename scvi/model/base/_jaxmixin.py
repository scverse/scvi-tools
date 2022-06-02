from typing import Optional, Union

import numpy as np

from scvi.dataloaders import DataSplitter
from scvi.train import JaxModuleInit, JaxTrainingPlan, TrainRunner


class JaxTrainingMixin:
    def train(
        self,
        max_epochs: Optional[int] = None,
        use_gpu: Optional[Union[str, int, bool]] = None,
        train_size: float = 0.9,
        validation_size: Optional[float] = None,
        batch_size: int = 128,
        lr: float = 1e-3,
        **trainer_kwargs,
    ):
        if max_epochs is None:
            n_cells = self.adata.n_obs
            max_epochs = np.min([round((20000 / n_cells) * 400), 400])

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
            self.module, use_gpu=use_gpu, optim_kwargs=dict(learning_rate=lr)
        )
        if "callbacks" not in trainer_kwargs.keys():
            trainer_kwargs["callbacks"] = []
        trainer_kwargs["callbacks"].append(JaxModuleInit())

        runner = TrainRunner(
            self,
            training_plan=self.training_plan,
            data_splitter=data_splitter,
            max_epochs=max_epochs,
            use_gpu=False,
            **trainer_kwargs,
        )
        runner()

        self.is_trained_ = True
        self.module.eval()
