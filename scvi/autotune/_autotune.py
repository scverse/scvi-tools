import logging
import os
import warnings
from typing import Callable, List, Optional, Union

import ray
import torch
from ray import tune
from ray.tune import ExperimentAnalysis
from ray.tune.schedulers import ASHAScheduler, PopulationBasedTraining
from ray.tune.schedulers.trial_scheduler import TrialScheduler

from scvi import model
from scvi._compat import Literal
from scvi._types import AnnOrMuData
from scvi.model.base import BaseModelClass

from ._utils import fetch_config, train_model

logger = logging.getLogger(__name__)


class ModelTuner:
    """
    Automated and parallel hyperparameter searches with Ray Tune. Wraps a :class:`~ray.tune.Tuner` instance that is
    attached to a scvi-tools model.

    Parameters
    ----------
    model_cls
        :class:`~scvi.model.base.BaseModelClass` on which to tune hyperparameters.
    adata
        AnnData/MuData object that has been registered via the model's ``setup_anndata`` or ``setup_mudata`` method.
    train_metrics:
        Metrics to compute during training.
    val_metrics:
        Metrics to compute after training.
    model_config:
        Configuration dictionary where keys correspond to keyword args in the model's ``__init__`` method.
        Values can be constants or instances of :class:`~ray.tune.search.sample.Domain`.
    trainer_config:
        Configuration dictionary where keys correspond to keyword args in the model's ``train`` method.
        Values can beconstants or instances of :class:`~ray.tune.search.sample.Domain`.
    plan_config:
        Configuration dictionary where keys correspond to keyword args in the model's :class:`~scvi.train.TrainingPlan`.
        Values can be constants or instances of :class:`~ray.tune.search.sample.Domain`.

    Notes
    -----
    See usage examples in the following tutorial(s):

    1. :doc:`/tutorials/notebooks/autotune`
    """

    def __init__(
        self,
        model_cls: BaseModelClass,
        adata: AnnOrMuData,
        train_metrics: Optional[List[Union[str, Callable]]] = None,
        val_metrics: Optional[List[Union[str, Callable]]] = None,
        model_config: Optional[dict] = None,
        trainer_config: Optional[dict] = None,
        plan_config: Optional[dict] = None,
    ):
        try:
            from ray import tune
        except ImportError:
            raise ImportError("Please install ray via `pip install ray`.")

        if model_cls not in [model.SCVI]:
            raise NotImplementedError(
                f"Tuning currently unsupported for {str(model_cls)}."
            )

        self.model_cls = model_cls
        self.adata = adata
        self.model_config = model_config or {}
        self.train_metrics = train_metrics or []
        self.val_metrics = val_metrics or []
        self.trainer_config = trainer_config or {}
        self.plan_config = plan_config or {}

        # Change fixed values to tune.choice with a single option
        for config in [self.model_config, self.trainer_config, self.plan_config]:
            for k, v in config:
                if not isinstance(v, tune.search.sample.Domain):
                    config[k] = tune.choice([v])

    def fit(
        self,
        metric,
        scheduler: Literal[
            "asha", "pbt", "pbt_replay", "pb_bandits", "hyperband", "bohb", "median"
        ] = "asha",
        resources_per_trial: Optional[dict] = None,
        **scheduler_kwargs,
    ) -> None:
        """
        Wrapper around `fit` of the current :class:`~ray.tune.Tuner` instance. Executes a hyperparameter search
        as specified.

        Parameters
        ----------
        metric
            Metric to optimize over. Must be present in ``self.train_metrics`` or ``self.val_metrics``.
        scheduler
            Ray Tune scheduler for trials. Options are:

            * ``'asha'``: `~ray.tune.schedulers.ASHAScheduler`.
            * ``'pbt'``: `~ray.tune.schedulers.PopulationBasedTraining`.
            * ``'pbt_replay'``: `~ray.tune.schedulers.PopulationBasedTrainingReplay`.
            * ``'pb_bandits'``: `~ray.tune.schedulers.pb2.PB2`.
            * ``'hyperband'``: `~ray.tune.schedulers.HyperBandScheduler`.
            * ``'bohb'``: `~ray.tune.schedulers.HyperBandForBOHB`.
            * ``'median'``: `~ray.tune.schedulers.MedianStoppingRule`.
        **scheduler_kwargs
            Keyword args for the scheduler. If not provided, appropriate defaults are used for each scheduler.
        """
        if metric not in self.train_metrics + self.val_metrics:
            raise ValueError(
                f"Metric {metric} not in tuner's train_metrics or val_metrics."
            )
        if scheduler == "asha":
            schedule = ASHAScheduler()
        elif scheduler == "pbt":
            schedule = PopulationBasedTraining()
        else:
            raise ValueError(
                f"Tuning currently unsupported for scheduler {scheduler}. Must be one of ['asha', 'pbt']."
            )
        print(schedule)

    def run(
        self,
        metric: str = None,
        scheduler: TrialScheduler = None,
        mode: str = "min",
        name: str = "scvi-experiment",
        num_samples: int = 1,
        resources_per_trial: Optional[dict] = None,
        local_dir: str = "./ray_results",
        train: bool = False,
        **kwargs,
    ) -> Union[BaseModelClass, ExperimentAnalysis]:
        """
        Wrap `tune.run`.

        Searches for the configuration of model, trainer, and training_plan
        hyperparameters that minimize or maximize the provided metric.
        Parameters
        ----------
        metric
            Metric to optimize over in self.metrics or from self.training_funcs
        scheduler
            Ray tune scheduler for trials. If `None`, defaults to
            :class:`~ray.tune.schedulers.ASHAScheduler`, with `max_t` set to the number
            of epochs.
        mode
            "min" or "max" to maximize or minimize the objective metric
        name
            Name of this experiment.
        num_samples
            Number of times to sample hyperparameters from the configuration space
        resources_per_trial
            Dictionary specifying the number of `gpu` and `cpu` in the optimization
        local_dir
            Local dir to save training results to.
        train
            Whether to train the resulting best model. Defaults to `False`
        **kwargs
            Additional arguments for func:`ray.tune.run`

        Returns
        -------
        A tuple with the best model object and tune Analysis object
        """
        if not scheduler:
            scheduler = ASHAScheduler(
                max_t=self.num_epochs, grace_period=1, reduction_factor=2
            )

        if not resources_per_trial:
            if torch.cuda.is_available():
                resources_per_trial = {"cpu": 1, "gpu": 1}
            else:
                resources_per_trial = {"cpu": 1}

        if resources_per_trial["cpu"] > 1 and resources_per_trial["gpu"]:
            resources_per_trial["gpu"] /= resources_per_trial["cpu"]

        if (
            metric in self.metrics
            and (
                "continuous_covariate_keys"
                or "categorical_covariate_keys" in self.setup_args
            )
            or self.top_hvg
        ):
            warnings.warn(
                f"You are optimizing over {metric} and testing different model architectures. "
                "This is not a recommended approach as the metric will be influenced"
                " by the architecture and not the model hyperparameters. Consider "
                "optimizing over a non-training metric, such as `autotune.metrics.silhouette_score`"
            )

        ray.init(num_cpus=resources_per_trial["cpu"], ignore_reinit_error=True)
        analysis = tune.run(
            self._trainable,
            metric=metric,
            mode=mode,
            config=self.config,
            num_samples=num_samples,
            scheduler=scheduler,
            progress_reporter=self.reporter,
            name=name,
            resources_per_trial=resources_per_trial,
            local_dir=local_dir,
            raise_on_failed_trial=False,
            **kwargs,
        )
        ray.shutdown()
        logger.info(
            "Best hyperparameters found were: {}".format(analysis.get_best_config())
        )

        best_config = analysis.best_config

        _, trainer_config, plan_config = fetch_config(self, best_config)

        # retrieve and load best model
        best_checkpoint = analysis.best_checkpoint
        print(best_checkpoint)
        best_model = self.model_cls.load(
            dir_path=os.path.join(best_checkpoint, "checkpoint"), adata=self.adata
        )
        best_model = train_model(
            self,
            best_model,
            trainer_config,
            plan_config,
            return_model=True,
            callbacks=False,
        )

        return best_model, analysis
