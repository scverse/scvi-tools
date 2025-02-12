from __future__ import annotations

import os
import warnings
from collections.abc import Callable
from copy import deepcopy
from datetime import datetime
from shutil import rmtree
from typing import TYPE_CHECKING

import flax
import numpy as np
import pyro
import torch
from anndata import AnnData
from lightning.pytorch.callbacks import Callback, ModelCheckpoint
from lightning.pytorch.callbacks.early_stopping import EarlyStopping
from lightning.pytorch.utilities import rank_zero_info

from scvi import settings
from scvi.model.base import BaseModelClass
from scvi.model.base._save_load import _load_saved_files
from scvi.utils import dependencies

if TYPE_CHECKING:
    from typing import Literal

    import lightning.pytorch as pl

    from scvi.dataloaders import AnnDataLoader

MetricCallable = Callable[[BaseModelClass], float]


class SaveCheckpoint(ModelCheckpoint):
    """``BETA`` Saves model checkpoints based on a monitored metric.

    Inherits from :class:`~lightning.pytorch.callbacks.ModelCheckpoint` and modifies the default
    behavior to save the full model state instead of just the state dict. This is necessary for
    compatibility with :class:`~scvi.model.base.BaseModelClass`.

    The best model save and best model score based on ``monitor`` can be accessed post-training
    with the ``best_model_path`` and ``best_model_score`` attributes, respectively.

    Known issues:

    * Does not set ``train_indices``, ``validation_indices``, and ``test_indices`` for checkpoints.
    * Does not set ``history`` for checkpoints. This can be accessed in the final model however.
    * Unsupported arguments: ``save_weights_only`` and ``save_last``.

    Parameters
    ----------
    dirpath
        Base directory to save the model checkpoints. If ``None``, defaults to a subdirectory in
        :attr:``scvi.settings.logging_dir`` formatted with the current date, time, and monitor.
    filename
        Name for the checkpoint directories, which can contain formatting options for auto-filling.
        If ``None``, defaults to ``{epoch}-{step}-{monitor}``.
    monitor
        Metric to monitor for checkpointing.
    load_best_on_end
        If ``True``, loads the best model state into the model at the end of training.
    **kwargs
        Additional keyword arguments passed into the constructor for
        :class:`~lightning.pytorch.callbacks.ModelCheckpoint`.
    """

    def __init__(
        self,
        dirpath: str | None = None,
        filename: str | None = None,
        monitor: str = "validation_loss",
        load_best_on_end: bool = False,
        **kwargs,
    ):
        if dirpath is None:
            dirpath = os.path.join(
                settings.logging_dir,
                datetime.now().strftime("%Y-%m-%d_%H-%M-%S") + f"_{monitor}",
            )
        if filename is None:
            filename = "{epoch}-{step}-{" + monitor + "}"
        if "save_weights_only" in kwargs:
            warnings.warn(
                "`save_weights_only` is not supported in `SaveCheckpoint` and will be ignored.",
                RuntimeWarning,
                stacklevel=settings.warnings_stacklevel,
            )
            kwargs.pop("save_weights_only")
        if "save_last" in kwargs:
            warnings.warn(
                "`save_last` is not supported in `SaveCheckpoint` and will be ignored.",
                RuntimeWarning,
                stacklevel=settings.warnings_stacklevel,
            )
            kwargs.pop("save_last")
        self.load_best_on_end = load_best_on_end

        super().__init__(
            dirpath=dirpath,
            filename=filename,
            monitor=monitor,
            **kwargs,
        )

    def on_save_checkpoint(self, trainer: pl.Trainer, *args) -> None:
        """Saves the model state on Lightning checkpoint saves."""
        # set post training state before saving
        model = trainer._model
        model.module.eval()
        model.is_trained_ = True
        model.trainer = trainer

        monitor_candidates = self._monitor_candidates(trainer)
        save_path = self.format_checkpoint_name(monitor_candidates)
        # by default, the function above gives a .ckpt extension
        save_path = save_path.split(".ckpt")[0]
        model.save(save_path, save_anndata=False, overwrite=True)

        model.module.train()
        model.is_trained_ = False

    def _remove_checkpoint(self, trainer: pl.Trainer, filepath: str) -> None:
        """Removes model saves that are no longer needed.

        Calls the superclass method and then removes the :class:`~scvi.model.base.BaseModelClass`
        save directory.
        """
        super()._remove_checkpoint(trainer, filepath)

        model_path = filepath.split(".ckpt")[0]
        if os.path.exists(model_path) and os.path.isdir(model_path):
            rmtree(model_path)

    def _update_best_and_save(
        self,
        current: torch.Tensor,
        trainer: pl.Trainer,
        monitor_candidates: dict[str, torch.Tensor],
    ) -> None:
        """Replaces Lightning checkpoints with :class:`~scvi.model.base.BaseModelClass` saves.

        Calls the superclass method and then replaces the Lightning checkpoint file with
        the :class:`~scvi.model.base.BaseModelClass` save directory.
        """
        super()._update_best_and_save(current, trainer, monitor_candidates)

        if os.path.exists(self.best_model_path):
            os.remove(self.best_model_path)
        self.best_model_path = self.best_model_path.split(".ckpt")[0]

    def on_train_end(self, trainer: pl.Trainer, pl_module: pl.LightningModule) -> None:
        """Loads the best model state into the model at the end of training."""
        if not self.load_best_on_end:
            return

        _, _, best_state_dict, _ = _load_saved_files(
            self.best_model_path,
            load_adata=False,
            map_location=pl_module.module.device,
        )
        pyro_param_store = best_state_dict.pop("pyro_param_store", None)
        pl_module.module.load_state_dict(best_state_dict)
        if pyro_param_store is not None:
            # For scArches shapes are changed and we don't want to overwrite these changed shapes.
            pyro.get_param_store().set_state(pyro_param_store)


class SubSampleLabels(Callback):
    """Subsample labels."""

    def __init__(self):
        super().__init__()

    def on_train_epoch_start(self, trainer, pl_module):
        """Subsample labels at the beginning of each epoch."""
        trainer.train_dataloader.resample_labels()
        super().on_train_epoch_start(trainer, pl_module)


class SaveBestState(Callback):
    r"""``DEPRECATED`` Save the best module state and restore into model.

    Parameters
    ----------
    monitor
        quantity to monitor.
    verbose
        verbosity, True or False.
    mode
        one of ["min", "max"].
    period
        Interval (number of epochs) between checkpoints.

    Examples
    --------
    from scvi.train import Trainer
    from scvi.train import SaveBestState

    Notes
    -----
    Lifecycle: deprecated in v1.2 and to be removed in v1.3. Please use
        :class:`~scvi.train.callbacks.SaveCheckpoint` instead.
    """

    def __init__(
        self,
        monitor: str = "elbo_validation",
        mode: str = "min",
        verbose=False,
        period=1,
    ):
        super().__init__()

        warnings.warn(
            "`SaveBestState` is deprecated in v1.2 and will be removed in v1.3. Please use "
            "`SaveCheckpoint` instead. See https://github.com/scverse/scvi-tools/issues/2568 "
            "for more details.",
            DeprecationWarning,
            stacklevel=settings.warnings_stacklevel,
        )

        self.monitor = monitor
        self.verbose = verbose
        self.period = period
        self.epochs_since_last_check = 0
        self.best_module_state = None

        if mode not in ["min", "max"]:
            raise ValueError(
                f"SaveBestState mode {mode} is unknown",
            )

        if mode == "min":
            self.monitor_op = np.less
            self.best_module_metric_val = np.inf
            self.mode = "min"
        elif mode == "max":
            self.monitor_op = np.greater
            self.best_module_metric_val = -np.inf
            self.mode = "max"
        else:
            if "acc" in self.monitor or self.monitor.startswith("fmeasure"):
                self.monitor_op = np.greater
                self.best_module_metric_val = -np.inf
                self.mode = "max"
            else:
                self.monitor_op = np.less
                self.best_module_metric_val = np.inf
                self.mode = "min"

    def check_monitor_top(self, current):
        return self.monitor_op(current, self.best_module_metric_val)

    def on_validation_epoch_end(self, trainer, pl_module):
        logs = trainer.callback_metrics
        self.epochs_since_last_check += 1

        if trainer.current_epoch > 0 and self.epochs_since_last_check >= self.period:
            self.epochs_since_last_check = 0
            current = logs.get(self.monitor)

            if current is None:
                warnings.warn(
                    f"Can save best module state only with {self.monitor} available, skipping.",
                    RuntimeWarning,
                    stacklevel=settings.warnings_stacklevel,
                )
            else:
                if isinstance(current, torch.Tensor):
                    current = current.item()
                if self.check_monitor_top(current):
                    self.best_module_state = deepcopy(pl_module.module.state_dict())
                    self.best_module_metric_val = current

                    if self.verbose:
                        rank_zero_info(
                            f"\nEpoch {trainer.current_epoch:05d}: {self.monitor} reached."
                            f" Module best state updated."
                        )

    def on_train_start(self, trainer, pl_module):
        self.best_module_state = deepcopy(pl_module.module.state_dict())

    def on_train_end(self, trainer, pl_module):
        pl_module.module.load_state_dict(self.best_module_state)


class LoudEarlyStopping(EarlyStopping):
    """Loud early stopping callback.

    Wrapper of :class:`~lightning.pytorch.callbacks.early_stopping.EarlyStopping callback that
    prints the reason for stopping on teardown. When the early stopping condition is met, the
    reason is saved to the callback instance, then printed on teardown. By printing on teardown, we
    do not interfere with the progress bar callback.
    """

    def __init__(self, **kwargs) -> None:
        super().__init__(**kwargs)
        self.early_stopping_reason = None

    def _evaluate_stopping_criteria(self, current: torch.Tensor) -> tuple[bool, str]:
        should_stop, reason = super()._evaluate_stopping_criteria(current)
        if should_stop:
            self.early_stopping_reason = reason
        return should_stop, reason

    def teardown(
        self,
        _trainer: pl.Trainer,
        _pl_module: pl.LightningModule,
        stage: str | None = None,
    ) -> None:
        """Print the reason for stopping on teardown."""
        if self.early_stopping_reason is not None:
            print(self.early_stopping_reason)


class JaxModuleInit(Callback):
    """A callback to initialize the Jax-based module."""

    def __init__(self, dataloader: AnnDataLoader = None) -> None:
        super().__init__()
        self.dataloader = dataloader

    def on_train_start(self, trainer, pl_module):
        module = pl_module.module
        if self.dataloader is None:
            dl = trainer.datamodule.train_dataloader()
        else:
            dl = self.dataloader
        module_init = module.init(module.rngs, next(iter(dl)))
        state, params = flax.core.pop(module_init, "params")
        pl_module.set_train_state(params, state)


@dependencies("scib_metrics")
class ScibCallback(Callback):
    # example to use in debug of the early stopping callback:
    # tune_callback = ScibCallback(stage="validation", metric="BioConservation")
    # stage="validation"
    # metric = "Bio conservation"
    # trainer=_trainer
    # pl_module=_pl_module
    from scib_metrics.benchmark import BatchCorrection, BioConservation

    def __init__(
        self,
        bio_conservation_metrics: BioConservation | None = None,
        batch_correction_metrics: BatchCorrection | None = None,
        stage: Literal["training", "validation", "both"] = "validation",
        metric: str | None = "Total",
        num_rows_to_select: int = 100,
    ):
        super().__init__()
        self.bio_conservation_metrics = bio_conservation_metrics
        self.batch_correction_metrics = batch_correction_metrics
        self.stage = stage
        self.metric = metric
        self.num_rows_to_select = num_rows_to_select

    def compute_metrics(
        self,
        trainer: pl.Trainer,
        pl_module: pl.LightningModule,
        stage: Literal["training", "validation"],
    ):
        from scib_metrics.benchmark import BatchCorrection, Benchmarker, BioConservation
        from scib_metrics.benchmark._core import metric_name_cleaner

        if self.metric is None:
            return
        if stage == "training" and self.stage not in ["training", "both"]:
            return
        elif stage == "validation" and self.stage not in ["validation", "both"]:
            return

        if not hasattr(pl_module, f"_{stage}_epoch_outputs"):
            raise ValueError(f"The training plan must have a `_{stage}_epoch_outputs` attribute.")

        outputs = getattr(pl_module, f"_{stage}_epoch_outputs")
        x = outputs["x"].numpy()
        # x = np.zeros(x.shape) #TODO: should we do it? can be done also in trainingplans already
        z = outputs["z"].numpy()
        batch = outputs["batch"].numpy()
        labels = outputs["labels"].numpy()

        # TODO: subsample to save time
        # rand_idx = np.random.choice(x.shape[0], self.num_rows_to_select, replace=False)
        # batch = batch[rand_idx]
        # labels = labels[rand_idx]
        # x = x[rand_idx]
        # z = z[rand_idx]

        # adjust which metric to run exactly
        found_metric = next(
            (key for key, value in metric_name_cleaner.items() if value == self.metric), None
        )
        # specal cases:
        if self.metric == "Leiden NMI" or self.metric == "Leiden ARI":
            found_metric = "nmi_ari_cluster_labels_leiden"
        if self.metric == "KMeans NMI" or self.metric == "KMeans ARI":
            found_metric = "nmi_ari_cluster_labels_kmeans"
        if found_metric is not None:
            # beucase originaly those classes are frozen we cant just set the metric to True
            # Need to do it manualy unfortunatley
            if found_metric == "isolated_labels":
                self.bio_conservation_metrics = BioConservation(True, False, False, False, False)
                self.batch_correction_metrics = BatchCorrection(False, False, False, False, False)
            if found_metric == "nmi_ari_cluster_labels_leiden":
                self.bio_conservation_metrics = BioConservation(False, True, False, False, False)
                self.batch_correction_metrics = BatchCorrection(False, False, False, False, False)
            if found_metric == "nmi_ari_cluster_labels_kmeans":
                self.bio_conservation_metrics = BioConservation(False, False, True, False, False)
                self.batch_correction_metrics = BatchCorrection(False, False, False, False, False)
            if found_metric == "silhouette_label":
                self.bio_conservation_metrics = BioConservation(False, False, False, True, False)
                self.batch_correction_metrics = BatchCorrection(False, False, False, False, False)
            if found_metric == "clisi_knn":
                self.bio_conservation_metrics = BioConservation(False, False, False, False, True)
                self.batch_correction_metrics = BatchCorrection(False, False, False, False, False)
            if found_metric == "silhouette_batch":
                self.bio_conservation_metrics = BioConservation(False, False, False, False, False)
                self.batch_correction_metrics = BatchCorrection(True, False, False, False, False)
            if found_metric == "ilisi_knn":
                self.bio_conservation_metrics = BioConservation(False, False, False, False, False)
                self.batch_correction_metrics = BatchCorrection(False, True, False, False, False)
            if found_metric == "kbet_per_label":
                self.bio_conservation_metrics = BioConservation(False, False, False, False, False)
                self.batch_correction_metrics = BatchCorrection(False, False, True, False, False)
            if found_metric == "graph_connectivity":
                self.bio_conservation_metrics = BioConservation(False, False, False, False, False)
                self.batch_correction_metrics = BatchCorrection(False, False, False, True, False)
            if found_metric == "pcr_comparison":
                self.bio_conservation_metrics = BioConservation(False, False, False, False, False)
                self.batch_correction_metrics = BatchCorrection(False, False, False, False, True)
        else:
            # its an aggregative metric
            if self.metric == "Total":
                # we jsut run them all, which is the default
                pass
            elif self.metric == "Batch correction":
                # we run all batch correction and no bio conservation
                self.bio_conservation_metrics = BioConservation(False, False, False, False, False)
            elif self.metric == "Bio conservation":
                # we run all bio conservation and no batch corredction
                self.batch_correction_metrics = BatchCorrection(False, False, False, False, False)
            else:
                # an invalid metric!
                raise ValueError(f"`{self.metric}` is an invalid metric in scib-metrics autotune.")

        adata = AnnData(X=x, obs={"batch": batch, "labels": labels}, obsm={"z": z})
        benchmarker = Benchmarker(
            adata,
            batch_key="batch",
            label_key="labels",
            embedding_obsm_keys=["z"],
            bio_conservation_metrics=self.bio_conservation_metrics,
            batch_correction_metrics=self.batch_correction_metrics,
        )
        benchmarker.benchmark()
        results = benchmarker.get_results(min_max_scale=False).to_dict()
        metrics = {f"training {self.metric}": results[self.metric]["z"]}
        pl_module.logger.log_metrics(metrics, trainer.global_step)

        delattr(pl_module, f"_{stage}_epoch_outputs")

    def on_train_epoch_end(self, trainer: pl.Trainer, pl_module: pl.LightningModule):
        self.compute_metrics(trainer, pl_module, "training")

    def on_validation_epoch_end(self, trainer: pl.Trainer, pl_module: pl.LightningModule):
        self.compute_metrics(trainer, pl_module, "validation")
