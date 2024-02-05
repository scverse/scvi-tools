from __future__ import annotations

from lightning.pytorch import LightningDataModule

from scvi._types import Tunable
from scvi.dataloaders import DataSplitter
from scvi.model._utils import get_max_epochs_heuristic, use_distributed_sampler
from scvi.train import TrainingPlan, TrainRunner
from scvi.utils._docstrings import devices_dsp


class UnsupervisedTrainingMixin:
    """General purpose unsupervised train method."""

    _data_splitter_cls = DataSplitter
    _training_plan_cls = TrainingPlan
    _train_runner_cls = TrainRunner

    @devices_dsp.dedent
    def train(
        self,
        max_epochs: int | None = None,
        accelerator: str = "auto",
        devices: int | list[int] | str = "auto",
        train_size: float = 0.9,
        validation_size: float | None = None,
        shuffle_set_split: bool = True,
        load_sparse_tensor: bool = False,
        batch_size: Tunable[int] = 128,
        early_stopping: bool = False,
        datasplitter_kwargs: dict | None = None,
        plan_kwargs: dict | None = None,
        data_module: LightningDataModule | None = None,
        **trainer_kwargs,
    ):
        """Train the model.

        Parameters
        ----------
        max_epochs
            Number of passes through the dataset. If ``None``, defaults to
            `np.min([round((20000 / n_cells) * 400), 400])`. Must be passed if ``data_module`` is
            passed in and does not have ``n_obs`` attribute.
        %(param_accelerator)s
        %(param_devices)s
        train_size
            Size of training set in the range ``[0.0, 1.0]``. Not used if ``data_module`` is passed
            in.
        validation_size
            Size of the test set. If ``None``, defaults to ``1 - train_size``. If
            ``train_size + validation_size < 1``, the remaining cells belong to a test set. Not used
            if ``data_module`` is passed in.
        shuffle_set_split
            Whether to shuffle indices before splitting. If ``False``, the val, train, and test set
            are split in the sequential order of the data according to ``validation_size`` and
            ``train_size`` percentages. Not used if ``data_module`` is passed in.
        load_sparse_tensor
            ``EXPERIMENTAL`` If ``True``, loads data with sparse CSR or CSC layout as a
            :class:`~torch.Tensor` with the same layout. Can lead to speedups in data transfers to
            GPUs, depending on the sparsity of the data. Not used if ``data_module`` is passed in.
        batch_size
            Minibatch size to use during training. Not used if ``data_module`` is passed in.
        early_stopping
            Perform early stopping. Additional arguments can be passed in `**kwargs`.
            See :class:`~scvi.train.Trainer` for further options.
        datasplitter_kwargs
            Additional keyword arguments passed into :class:`~scvi.dataloaders.DataSplitter`. Not
            used if ``data_module`` is passed in.
        plan_kwargs
            Keyword args for :class:`~scvi.train.TrainingPlan`. Keyword arguments passed to
            `train()` will overwrite values present in `plan_kwargs`, when appropriate.
        data_module
            ``EXPERIMENTAL`` A :class:`~lightning.pytorch.LightningDataModule` instance to use for
            training. Can only be passed in if the model was not initialized with
            :class:`~anndata.AnnData`.
        **trainer_kwargs
            Other keyword args for :class:`~scvi.train.Trainer`.
        """
        if data_module is not None and not self._module_init_on_train:
            raise ValueError(
                "Cannot pass in `data_module` if the model was initialized with `adata`."
            )
        elif data_module is None and self._module_init_on_train:
            raise ValueError(
                "If the model was not initialized with `adata`, a `data_module` must be passed in."
            )

        if max_epochs is None and data_module is None:
            max_epochs = get_max_epochs_heuristic(self.adata.n_obs)
        elif (
            max_epochs is None
            and data_module is not None
            and hasattr(data_module, "n_obs")
        ):
            max_epochs = get_max_epochs_heuristic(data_module.n_obs)
        elif max_epochs is None and data_module is not None:
            raise ValueError(
                "If `data_module` does not have `n_obs` attribute, `max_epochs` must be passed in."
            )

        plan_kwargs = plan_kwargs or {}
        datasplitter_kwargs = datasplitter_kwargs or {}

        if data_module is None:
            data_module = self._data_splitter_cls(
                self.adata_manager,
                train_size=train_size,
                validation_size=validation_size,
                batch_size=batch_size,
                shuffle_set_split=shuffle_set_split,
                distributed_sampler=use_distributed_sampler(
                    trainer_kwargs.get("strategy", None)
                ),
                load_sparse_tensor=load_sparse_tensor,
                **datasplitter_kwargs,
            )
        else:
            self.module = self._module_cls(
                data_module.n_vars,
                n_batch=data_module.n_batch,
                n_labels=getattr(data_module, "n_labels", 1),
                n_continuous_cov=getattr(data_module, "n_continuous_cov", 0),
                n_cats_per_cov=getattr(data_module, "n_cats_per_cov", None),
                **self._module_kwargs,
            )
        training_plan = self._training_plan_cls(self.module, **plan_kwargs)

        es = "early_stopping"
        trainer_kwargs[es] = (
            early_stopping if es not in trainer_kwargs.keys() else trainer_kwargs[es]
        )
        runner = self._train_runner_cls(
            self,
            training_plan=training_plan,
            data_splitter=data_module,
            max_epochs=max_epochs,
            accelerator=accelerator,
            devices=devices,
            **trainer_kwargs,
        )
        return runner()
