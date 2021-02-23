import logging
from typing import Optional, Union

import pytorch_lightning as pl

from scvi.dataloaders import DataSplitter, SemiSupervisedDataSplitter
from scvi.lightning import Trainer
from scvi.model._utils import parse_use_gpu_arg
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
    use_gpu
        Use default GPU if available (if None or True), or index of GPU to use (if int),
        or name of GPU (if str), or use CPU (if False).
    trainer_kwargs
        Extra kwargs for :class:`~scvi.lightning.Trainer`

    Examples
    --------
    >>> # Following code should be within a subclass of BaseModelClass
    >>> data_splitter = DataSplitter(self.adata)
    >>> training_plan = TrainingPlan(self.module, len(data_splitter.train_idx))
    >>> runner = TrainRunner(
    >>>     self,
    >>>     training_plan=trianing_plan,
    >>>     data_splitter=data_splitter,
    >>>     max_epochs=max_epochs)
    >>> runner()
    """

    def __init__(
        self,
        model: BaseModelClass,
        training_plan: pl.LightningModule,
        data_splitter: Union[SemiSupervisedDataSplitter, DataSplitter],
        max_epochs: int,
        use_gpu: Optional[Union[str, int, bool]] = None,
        **trainer_kwargs,
    ):
        self.training_plan = training_plan
        self.data_splitter = data_splitter
        self.model = model
        gpus, device = parse_use_gpu_arg(use_gpu)
        self.gpus = gpus
        self.device = device
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
        self.model.is_trained_ = True
        self.model.to_device(self.device)
        self.model.trainer = self.trainer
