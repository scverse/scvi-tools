from typing import Optional

import numpy as np
import torch
from anndata import AnnData
from sklearn.model_selection._split import _validate_shuffle_split

from scvi import settings
from scvi.dataloaders._ann_dataloader import AnnDataLoader
from scvi.lightning import Trainer, TrainingPlan


class TrainRunner:
    def __init__(self, model_class, trainer, training_plan, data_splitter):
        self.training_plan = training_plan
        self.data_splitter = data_splitter
        self.model_class = model_class
        self.trainer = trainer

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

        # if use_gpu:
        #     self.module.cuda()

        # self.is_trained_ = True


class UnsupervisedDataSplitter:
    def __init__(
        self,
        adata: AnnData,
        train_size: float = 0.9,
        validation_size: Optional[float] = None,
        **kwargs,
    ):
        self.adata = adata
        self.train_size = train_size
        self.validation_size = validation_size
        self.data_loader_kwargs = kwargs
        self.train_idx, self.test_idx, self.val_idx = self.make_splits()

    def make_splits(self):
        train_size = float(self.train_size)
        if train_size > 1.0 or train_size <= 0.0:
            raise ValueError(
                "train_size needs to be greater than 0 and less than or equal to 1"
            )
        try:
            n = self.adata.n_obs
            n_train, n_val = _validate_shuffle_split(
                n, self.validation_size, self.train_size
            )
        except ValueError:
            if train_size != 1.0:
                raise ValueError(
                    "Choice of train_size={} and validation_size={} not understood".format(
                        self.train_size, self.validation_size
                    )
                )
            n_train, n_val = n, 0
        random_state = np.random.RandomState(seed=settings.seed)
        permutation = random_state.permutation(n)
        val_idx = permutation[:n_val]
        train_idx = permutation[n_val : (n_val + n_train)]
        test_idx = permutation[(n_val + n_train) :]
        return train_idx, test_idx, val_idx

    def __call__(self, remake_splits=False):
        if remake_splits:
            self.make_splits

        # do not remove drop_last=3, skips over small minibatches
        return (
            AnnDataLoader(
                self.adata,
                indices=self.train_idx,
                shuffle=True,
                drop_last=3,
                **self.data_loader_kwargs,
            ),
            AnnDataLoader(
                self.adata,
                indices=self.val_idx,
                shuffle=True,
                drop_last=3,
                **self.data_loader_kwargs,
            ),
            AnnDataLoader(
                self.adata,
                indices=self.test_idx,
                shuffle=True,
                drop_last=3,
                **self.data_loader_kwargs,
            ),
        )


class SemisupervisedDataSplitter:
    def __init__(self):
        pass


class UnsupervisedTrainingMixin:
    def train(
        self,
        max_epochs: Optional[int] = None,
        use_gpu: Optional[bool] = None,
        train_size: float = 0.9,
        validation_size: Optional[float] = None,
        batch_size: int = 128,
        plan_kwargs: Optional[dict] = None,
        **trainer_kwargs,
    ):
        """
        Train the model.

        Parameters
        ----------
        max_epochs
            Number of passes through the dataset. If `None`, defaults to
            `np.min([round((20000 / n_cells) * 400), 400])`
        use_gpu
            If `True`, use the GPU if available. Will override the use_gpu option when initializing model
        train_size
            Size of training set in the range [0.0, 1.0].
        validation_size
            Size of the test set. If `None`, defaults to 1 - `train_size`. If
            `train_size + validation_size < 1`, the remaining cells belong to a test set.
        batch_size
            Minibatch size to use during training.
        plan_kwargs
            Keyword args for model-specific Pytorch Lightning task. Keyword arguments passed to
            `train()` will overwrite values present in `plan_kwargs`, when appropriate.
        plan_class
            Optional override to use a specific TrainingPlan-type class.
        **trainer_kwargs
            Other keyword args for :class:`~scvi.lightning.Trainer`.
        """
        if max_epochs is None:
            n_cells = self.adata.n_obs
            max_epochs = np.min([round((20000 / n_cells) * 400), 400])

        if use_gpu is None:
            use_gpu = torch.cuda.is_available()
        else:
            use_gpu = use_gpu and torch.cuda.is_available()
        gpus = 1 if use_gpu else None

        pin_memory = (
            True if (settings.dl_pin_memory_gpu_training and use_gpu) else False
        )
        plan_kwargs = plan_kwargs if isinstance(plan_kwargs, dict) else dict()

        data_splitter = UnsupervisedDataSplitter(
            self.adata,
            train_size=train_size,
            validation_size=validation_size,
            pin_memory=pin_memory,
            batch_size=batch_size,
        )
        training_plan = TrainingPlan(
            self.module, len(data_splitter.train_idx), **plan_kwargs
        )
        trainer = Trainer(max_epochs=max_epochs, gpus=gpus, **trainer_kwargs)
        runner = TrainRunner(self, trainer, training_plan, data_splitter)
        return runner()
