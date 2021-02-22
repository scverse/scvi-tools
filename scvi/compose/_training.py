import logging

from scvi.lightning import Trainer

logger = logging.getLogger(__name__)


class TrainRunner:
    def __init__(
        self,
        model_class,
        training_plan,
        data_splitter,
        max_epochs,
        gpus,
        **trainer_kwargs,
    ):
        self.training_plan = training_plan
        self.data_splitter = data_splitter
        self.model_class = model_class
        self.gpus = gpus
        self.trainer = Trainer(max_epochs=max_epochs, gpus=gpus, **trainer_kwargs)

    def __call__(self):
        train_dl, val_dl, test_dl = self.data_splitter()
        self.model_class.train_indices = train_dl.indices
        self.model_class.test_indices = test_dl.indices
        self.model_class.validation_indices = val_dl.indices

        if len(val_dl.indices) == 0:
            # circumvent the empty data loader problem if all dataset used for training
            self.trainer.fit(self.training_plan, train_dl)
        else:
            self.trainer.fit(self.training_plan, train_dl, val_dl)
        try:
            self.model_class.history_ = self.trainer.logger.history
        except AttributeError:
            self.history_ = None

        self.model_class.module.eval()

        if self.gpus != 0:
            self.model_class.module.cuda()

        self.model_class.is_trained_ = True
