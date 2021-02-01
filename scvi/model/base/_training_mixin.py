from typing import Optional
from anndata import AnnData

import logging
import numpy as np
import torch

from typing import Union
from scvi import settings
from scvi.lightning._callbacks import SubSampleLabels
from scvi.lightning import Trainer, AdversarialTrainingPlan
from sklearn.model_selection._split import _validate_shuffle_split
from scvi.dataloaders import SemiSupervisedDataLoader, AnnDataLoader
from scvi.lightning._trainingplans import SemiSupervisedTrainingPlan, TrainingPlan

logger = logging.getLogger(__name__)


class UnsupervisedTrainingMixin:
    def train(
        self,
        max_epochs: Optional[int] = None,
        use_gpu: Optional[bool] = None,
        train_size: float = 0.9,
        validation_size: Optional[float] = None,
        batch_size: int = 128,
        plan_kwargs: Optional[dict] = None,
        plan_class=TrainingPlan,
        **kwargs,
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
        **kwargs
            Other keyword args for :class:`~scvi.lightning.Trainer`.
        """
        if use_gpu is None:
            use_gpu = self.use_gpu
        else:
            use_gpu = use_gpu and torch.cuda.is_available()

        gpus = 1 if use_gpu else None
        pin_memory = (
            True if (settings.dl_pin_memory_gpu_training and use_gpu) else False
        )

        if max_epochs is None:
            n_cells = self.adata.n_obs
            max_epochs = np.min([round((20000 / n_cells) * 400), 400])

        self.trainer = Trainer(
            max_epochs=max_epochs,
            gpus=gpus,
            **kwargs,
        )
        train_dl, val_dl, test_dl = self._train_test_val_split(
            self.adata,
            train_size=train_size,
            validation_size=validation_size,
            pin_memory=pin_memory,
            batch_size=batch_size,
        )
        self.train_indices_ = train_dl.indices
        self.test_indices_ = test_dl.indices
        self.validation_indices_ = val_dl.indices

        plan_kwargs = plan_kwargs if isinstance(plan_kwargs, dict) else dict()
        self._pl_task = plan_class(
            vae_model=self.module,
            n_obs_training=len(self.train_indices_),
            **plan_kwargs,
        )

        if train_size == 1.0:
            # circumvent the empty data loader problem if all dataset used for training
            self.trainer.fit(self._pl_task, train_dl)
        else:
            self.trainer.fit(self._pl_task, train_dl, val_dl)
        try:
            self.history_ = self.trainer.logger.history
        except AttributeError:
            self.history_ = None
        self.module.eval()
        if use_gpu:
            self.module.cuda()
        self.is_trained_ = True

    def _train_test_val_split(
        self,
        adata: AnnData,
        train_size: float = 0.9,
        validation_size: Optional[float] = None,
        **kwargs,
    ):
        """
        Creates data loaders ``train_set``, ``validation_set``, ``test_set``.

        If ``train_size + validation_set < 1`` then ``test_set`` is non-empty.

        Parameters
        ----------
        adata
            Setup AnnData to be split into train, test, validation sets
        train_size
            float, or None (default is 0.9)
        validation_size
            float, or None (default is None)
        **kwargs
            Keyword args for `_make_scvi_dl()`
        """
        train_size = float(train_size)
        if train_size > 1.0 or train_size <= 0.0:
            raise ValueError(
                "train_size needs to be greater than 0 and less than or equal to 1"
            )

        n = len(adata)
        try:
            n_train, n_val = _validate_shuffle_split(n, validation_size, train_size)
        except ValueError:
            if train_size != 1.0:
                raise ValueError(
                    "Choice of train_size={} and validation_size={} not understood".format(
                        train_size, validation_size
                    )
                )
            n_train, n_val = n, 0
        random_state = np.random.RandomState(seed=settings.seed)
        permutation = random_state.permutation(n)
        indices_validation = permutation[:n_val]
        indices_train = permutation[n_val : (n_val + n_train)]
        indices_test = permutation[(n_val + n_train) :]

        return (
            self._make_scvi_dl(
                adata,
                indices=indices_train,
                shuffle=True,
                scvi_dl_class=AnnDataLoader,
                **kwargs,
            ),
            self._make_scvi_dl(
                adata,
                indices=indices_validation,
                shuffle=True,
                scvi_dl_class=AnnDataLoader,
                **kwargs,
            ),
            self._make_scvi_dl(
                adata,
                indices=indices_test,
                shuffle=True,
                scvi_dl_class=AnnDataLoader,
                **kwargs,
            ),
        )


class AdversarialTrainingMixin(UnsupervisedTrainingMixin):
    def train(
        self,
        max_epochs: Optional[int] = 400,
        lr: float = 4e-3,
        use_gpu: Optional[bool] = None,
        train_size: float = 0.9,
        validation_size: Optional[float] = None,
        batch_size: int = 256,
        early_stopping: bool = True,
        check_val_every_n_epoch: Optional[int] = None,
        reduce_lr_on_plateau: bool = True,
        n_steps_kl_warmup: Union[int, None] = None,
        n_epochs_kl_warmup: Union[int, None] = None,
        adversarial_classifier: Optional[bool] = None,
        plan_kwargs: Optional[dict] = None,
        **kwargs,
    ):
        """
        Trains the model using amortized variational inference.

        Parameters
        ----------
        max_epochs
            Number of passes through the dataset.
        lr
            Learning rate for optimization.
        use_gpu
            If `True`, use the GPU if available.
        train_size
            Size of training set in the range [0.0, 1.0].
        validation_size
            Size of the test set. If `None`, defaults to 1 - `train_size`. If
            `train_size + validation_size < 1`, the remaining cells belong to a test set.
        batch_size
            Minibatch size to use during training.
        early_stopping
            Whether to perform early stopping with respect to the validation set.
        check_val_every_n_epoch
            Check val every n train epochs. By default, val is not checked, unless `early_stopping` is `True`
            or `reduce_lr_on_plateau` is `True`. If either of the latter conditions are met, val is checked
            every epoch.
        reduce_lr_on_plateau
            Reduce learning rate on plateau of validation metric (default is ELBO).
        n_steps_kl_warmup
            Number of training steps (minibatches) to scale weight on KL divergences from 0 to 1.
            Only activated when `n_epochs_kl_warmup` is set to None. If `None`, defaults
            to `floor(0.75 * adata.n_obs)`.
        n_epochs_kl_warmup
            Number of epochs to scale weight on KL divergences from 0 to 1.
            Overrides `n_steps_kl_warmup` when both are not `None`.
        adversarial_classifier
            Whether to use adversarial classifier in the latent space. This helps mixing when
            there are missing proteins in any of the batches. Defaults to `True` is missing proteins
            are detected.
        plan_kwargs
            Keyword args for :class:`~scvi.lightning.AdversarialTrainingPlan`. Keyword arguments passed to
            `train()` will overwrite values present in `plan_kwargs`, when appropriate.
        **kwargs
            Other keyword args for :class:`~scvi.lightning.Trainer`.
        """
        if adversarial_classifier is None:
            imputation = (
                True if "totalvi_batch_mask" in self.scvi_setup_dict_.keys() else False
            )
            adversarial_classifier = True if imputation else False

        n_steps_kl_warmup = (
            n_steps_kl_warmup
            if n_steps_kl_warmup is not None
            else int(0.75 * self.adata.n_obs)
        )
        if reduce_lr_on_plateau:
            check_val_every_n_epoch = 1

        update_dict = {
            "lr": lr,
            "adversarial_classifier": adversarial_classifier,
            "reduce_lr_on_plateau": reduce_lr_on_plateau,
            "n_epochs_kl_warmup": n_epochs_kl_warmup,
            "n_steps_kl_warmup": n_steps_kl_warmup,
            "check_val_every_n_epoch": check_val_every_n_epoch,
        }
        if plan_kwargs is not None:
            plan_kwargs.update(update_dict)
        else:
            plan_kwargs = update_dict

        super().train(
            max_epochs=max_epochs,
            use_gpu=use_gpu,
            train_size=train_size,
            validation_size=validation_size,
            batch_size=batch_size,
            early_stopping=early_stopping,
            plan_kwargs=plan_kwargs,
            plan_class=AdversarialTrainingPlan,
            **kwargs,
        )


class SemiSupervisedTrainingMixin:
    def train(
        self,
        max_epochs: Optional[int] = None,
        n_samples_per_label: Optional[float] = None,
        check_val_every_n_epoch: Optional[int] = None,
        train_size: float = 0.9,
        validation_size: Optional[float] = None,
        batch_size: int = 128,
        use_gpu: Optional[bool] = None,
        plan_kwargs: Optional[dict] = None,
        **kwargs,
    ):
        """
        Train the model.

        Parameters
        ----------
        max_epochs
            Number of passes through the dataset for semisupervised training.
        n_samples_per_label
            Number of subsamples for each label class to sample per epoch. By default, there
            is no label subsampling.
        check_val_every_n_epoch
            Frequency with which metrics are computed on the data for validation set for both
            the unsupervised and semisupervised trainers. If you'd like a different frequency for
            the semisupervised trainer, set check_val_every_n_epoch in semisupervised_train_kwargs.
        train_size
            Size of training set in the range [0.0, 1.0].
        validation_size
            Size of the test set. If `None`, defaults to 1 - `train_size`. If
            `train_size + validation_size < 1`, the remaining cells belong to a test set.
        batch_size
            Minibatch size to use during training.
        use_gpu
            If `True`, use the GPU if available. Will override the use_gpu option when initializing model
        plan_kwargs
            Keyword args for :class:`~scvi.lightning.SemiSupervisedTrainingPlan`. Keyword arguments passed to
            `train()` will overwrite values present in `plan_kwargs`, when appropriate.
        **kwargs
            Other keyword args for :class:`~scvi.lightning.Trainer`.
        """
        if max_epochs is None:
            n_cells = self.adata.n_obs
            max_epochs = np.min([round((20000 / n_cells) * 400), 400])

        logger.info("Training for {} epochs.".format(max_epochs))

        use_gpu = use_gpu if use_gpu is not None else self.use_gpu
        gpus = 1 if use_gpu else None
        pin_memory = (
            True if (settings.dl_pin_memory_gpu_training and use_gpu) else False
        )
        train_dl, val_dl, test_dl = self._train_test_val_split(
            self.adata,
            unlabeled_category=self.unlabeled_category_,
            train_size=train_size,
            validation_size=validation_size,
            n_samples_per_label=n_samples_per_label,
            pin_memory=pin_memory,
            batch_size=batch_size,
        )

        self.train_indices_ = train_dl.indices
        self.validation_indices_ = val_dl.indices
        self.test_indices_ = test_dl.indices

        plan_kwargs = {} if plan_kwargs is None else plan_kwargs
        self._task = SemiSupervisedTrainingPlan(self.module, **plan_kwargs)

        # if we have labeled cells, we want to subsample labels each epoch
        sampler_callback = (
            [SubSampleLabels()] if len(self._labeled_indices) != 0 else []
        )
        self._trainer = Trainer(
            max_epochs=max_epochs,
            gpus=gpus,
            callbacks=sampler_callback,
            check_val_every_n_epoch=check_val_every_n_epoch,
            **kwargs,
        )
        if len(self.validation_indices_) != 0:
            self._trainer.fit(self._task, train_dl, val_dl)
        else:
            self._trainer.fit(self._task, train_dl)
        self.module.eval()
        self.is_trained_ = True

    def _train_test_val_split(
        self,
        adata: AnnData,
        unlabeled_category,
        train_size: float = 0.9,
        validation_size: Optional[float] = None,
        n_samples_per_label: Optional[int] = None,
        **kwargs,
    ):
        """
        Creates data loaders ``train_set``, ``validation_set``, ``test_set``.

        If ``train_size + validation_set < 1`` then ``test_set`` is non-empty.
        The ratio between labeled and unlabeled data in adata will be preserved
        in the train/test/val sets.

        Parameters
        ----------
        adata
            AnnData to split into train/test/val sets
        unlabeled_category
            Category to treat as unlabeled
        train_size
            float, or None (default is 0.9)
        validation_size
            float, or None (default is None)
        n_samples_per_label
            Number of subsamples for each label class to sample per epoch
        **kwargs
            Keyword args for `_make_scvi_dl()`
        """
        train_size = float(train_size)
        if train_size > 1.0 or train_size <= 0.0:
            raise ValueError(
                "train_size needs to be greater than 0 and less than or equal to 1"
            )

        n_labeled_idx = len(self._labeled_indices)
        n_unlabeled_idx = len(self._unlabeled_indices)

        def get_train_val_split(n_samples, test_size, train_size):
            try:
                n_train, n_val = _validate_shuffle_split(
                    n_samples, test_size, train_size
                )
            except ValueError:
                if train_size != 1.0 and n_samples != 1:
                    raise ValueError(
                        "Choice of train_size={} and validation_size={} not understood".format(
                            train_size, test_size
                        )
                    )
                n_train, n_val = n_samples, 0
            return n_train, n_val

        if n_labeled_idx != 0:
            n_labeled_train, n_labeled_val = get_train_val_split(
                n_labeled_idx, validation_size, train_size
            )
            labeled_permutation = np.random.choice(
                self._labeled_indices, len(self._labeled_indices), replace=False
            )
            labeled_idx_val = labeled_permutation[:n_labeled_val]
            labeled_idx_train = labeled_permutation[
                n_labeled_val : (n_labeled_val + n_labeled_train)
            ]
            labeled_idx_test = labeled_permutation[(n_labeled_val + n_labeled_train) :]
        else:
            labeled_idx_test = []
            labeled_idx_train = []
            labeled_idx_val = []

        if n_unlabeled_idx != 0:
            n_unlabeled_train, n_unlabeled_val = get_train_val_split(
                n_unlabeled_idx, validation_size, train_size
            )
            unlabeled_permutation = np.random.choice(
                self._unlabeled_indices, len(self._unlabeled_indices)
            )
            unlabeled_idx_val = unlabeled_permutation[:n_unlabeled_val]
            unlabeled_idx_train = unlabeled_permutation[
                n_unlabeled_val : (n_unlabeled_val + n_unlabeled_train)
            ]
            unlabeled_idx_test = unlabeled_permutation[
                (n_unlabeled_val + n_unlabeled_train) :
            ]
        else:
            unlabeled_idx_train = []
            unlabeled_idx_val = []
            unlabeled_idx_test = []

        indices_train = np.concatenate((labeled_idx_train, unlabeled_idx_train))
        indices_val = np.concatenate((labeled_idx_val, unlabeled_idx_val))
        indices_test = np.concatenate((labeled_idx_test, unlabeled_idx_test))

        indices_train = indices_train.astype(int)
        indices_val = indices_val.astype(int)
        indices_test = indices_test.astype(int)

        if len(self._labeled_indices) != 0:
            dataloader_class = SemiSupervisedDataLoader
            dl_kwargs = {
                "unlabeled_category": unlabeled_category,
                "n_samples_per_label": n_samples_per_label,
            }
        else:
            dataloader_class = AnnDataLoader
            dl_kwargs = {}
        dl_kwargs.update(kwargs)

        scanvi_train_dl = self._make_scvi_dl(
            adata,
            indices=indices_train,
            shuffle=True,
            scvi_dl_class=dataloader_class,
            **dl_kwargs,
        )
        scanvi_val_dl = self._make_scvi_dl(
            adata,
            indices=indices_val,
            shuffle=True,
            scvi_dl_class=dataloader_class,
            **dl_kwargs,
        )
        scanvi_test_dl = self._make_scvi_dl(
            adata,
            indices=indices_test,
            shuffle=True,
            scvi_dl_class=dataloader_class,
            **dl_kwargs,
        )

        return scanvi_train_dl, scanvi_val_dl, scanvi_test_dl
