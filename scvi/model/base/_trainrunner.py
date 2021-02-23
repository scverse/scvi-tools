import logging
from typing import List, Union

import pytorch_lightning as pl

from scvi.dataloaders import DataSplitter, SemiSupervisedDataSplitter
from scvi.lightning import Trainer
from scvi.model.base import BaseModelClass

logger = logging.getLogger(__name__)


class TrainRunner:
    """
    TrainRunner calls Trainer.fit() and handles pre and post training procedures.

    Parameters
    ----------
    model
        model to train
    training_plan
        initialized TrainingPlan
    data_splitter
        initialized :class:`~scvi.dataloaders.SemiSupervisedDataSplitter` or
        :class:`~scvi.dataloaders.DataSplitter`
    max_epochs
        max_epochs to train for
    gpus
        Which GPU to train on
    trainer_kwargs
        Extra kwargs for :class:`~scvi.lightning.Trainer`
    """

    def __init__(
        self,
        model: BaseModelClass,
        training_plan: pl.LightningModule,
        data_splitter: Union[SemiSupervisedDataSplitter, DataSplitter],
        max_epochs: int,
        gpus: Union[List[int], str, int],
        **trainer_kwargs,
    ):
        self.training_plan = training_plan
        self.data_splitter = data_splitter
        self.model = model
        self.gpus = gpus
        self.trainer = Trainer(max_epochs=max_epochs, gpus=gpus, **trainer_kwargs)

    def __call__(self):
        train_dl, val_dl, test_dl = self.data_splitter()
        self.model.train_indices = train_dl.indices
        self.model.test_indices = test_dl.indices
        self.model.validation_indices = val_dl.indices

        if len(val_dl.indices) == 0:
            # circumvent the empty data loader problem if all dataset used for training
            self.trainer.fit(self.training_plan, train_dl)
        else:
            self.trainer.fit(self.training_plan, train_dl, val_dl)
        try:
            self.model.history_ = self.trainer.logger.history
        except AttributeError:
            self.history_ = None

        self.model.module.eval()

        if self.gpus != 0:
            self.model.module.cuda()

        self.model.is_trained_ = True
        self.model.trainer = self.trainer
