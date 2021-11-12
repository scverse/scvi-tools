from typing import Optional

import anndata
import torch
from ray import tune
from ray.tune import CLIReporter
from ray.tune.schedulers import ASHAScheduler
from ray.tune.schedulers.trial_scheduler import TrialScheduler
from ray.tune.suggest.hyperopt import HyperOptSearch
from ray.tune.suggest.search import SearchAlgorithm

from ._callbacks import ModelSave, _TuneReportMetricFunctionsCallback


class Autotune:
    """
    Hyperparameter tuning using Ray Tune.

    Parameters
    ----------
    adata
        AnnData object we will tune the model on.
    model
        Model from scvi.model we will tune.
    training_metrics
        Metrics to track during training.
    metric_functions
        For metrics calculated after training a model, like silhouette distance.
    model_hyperparams
        Config for the model hyperparameters https://docs.ray.io/en/master/tune/api_docs/search_space.html.
    trainer_hyperparams
        Config for the trainer hyperparameters https://docs.ray.io/en/master/tune/api_docs/search_space.html.
    plan_hyperparams
        Config for the training_plan hyperparameters https://docs.ray.io/en/master/tune/api_docs/search_space.html.
    """

    def __init__(
        self,
        adata: anndata.AnnData,
        model,
        training_metrics: Optional[list] = None,
        metric_functions: Optional[dict] = None,
        model_hyperparams: Optional[dict] = None,
        trainer_hyperparams: Optional[dict] = None,
        plan_hyperparams: Optional[dict] = None,
        num_epochs: int = 2,
    ):
        if not training_metrics:
            training_metrics = []
        if not metric_functions:
            metric_functions = {}
        if not model_hyperparams:
            model_hyperparams = {}
        if not trainer_hyperparams:
            trainer_hyperparams = {}
        if not plan_hyperparams:
            plan_hyperparams = {}
        self.adata = adata
        self.model = model
        self.training_metrics = training_metrics
        self.metric_functions = metric_functions
        self.model_hyperparams = model_hyperparams
        self.trainer_hyperparams = trainer_hyperparams
        self.plan_hyperparams = plan_hyperparams
        self.metrics = training_metrics
        self.reporter = CLIReporter(
            metric_columns=training_metrics + list(self.metric_functions.keys())
        )
        self.config = {}
        for d in [model_hyperparams, trainer_hyperparams, plan_hyperparams]:
            if d is not None:
                self.config.update(d)
        self.num_epochs = num_epochs

    def _trainable(self, config, checkpoint_dir=None):
        model_config = {}
        trainer_config = {}
        plan_config = {}
        for key in config:
            if key in self.model_hyperparams:
                model_config[key] = config[key]
            elif key in self.trainer_hyperparams:
                trainer_config[key] = config[key]
            elif key in self.plan_hyperparams:
                plan_config[key] = config[key]

        _model = self.model(self.adata, **model_config)
        _model.train(
            **trainer_config,
            plan_kwargs=plan_config,
            callbacks=[
                ModelSave(_model),
                _TuneReportMetricFunctionsCallback(
                    metrics=self.metrics,
                    on="validation_end",
                    model=_model,
                    metric_functions=self.metric_functions,
                ),
            ],
            check_val_every_n_epoch=1,
            max_epochs=self.num_epochs,
        )

    def run(
        self,
        metric: str = None,
        scheduler: TrialScheduler = None,
        search_alg: SearchAlgorithm = None,
        mode: str = "min",
        name: str = "scvi-experiment",
        num_samples: int = 10,
        resources_per_trial: Optional[dict] = None,
        **kwargs,
    ):
        """
        Wrapper for `tune.run`. Searches for the configuration of model, trainer, and training_plan
        hyperparameters that minimize or maximize the provided metric.

        Parameters
        ----------
        metric
            Metric to optimize over in self.metrics or from self.training_funcs
        scheduler
            Ray tune scheduler for trials defaults to ASHA.
        search_alg
            Search algorithm from `tune.suggest`, defaults to hyperopt.
        mode
            "min" or "max" to maximize or minimize the objective metric
        name
            Name of this experiment.
        num_samples
            Number of times to sample hyperparameters from the configuration space.
        **kwargs
            Additional arguments for `tune.run`

        Returns
        -------
        A tuple with the best model object and tune Analysis object

        """
        if not scheduler:
            scheduler = ASHAScheduler(max_t=2, grace_period=1, reduction_factor=2)

        if not search_alg:
            search_alg = HyperOptSearch(space=self.config, mode=mode, metric=metric)

        if not resources_per_trial:
            if torch.cuda.is_available():
                resources_per_trial = {"cpu": 1, "gpu": 1}
            else:
                resources_per_trial = {"cpu": 1}

        analysis = tune.run(
            self._trainable,
            metric=metric,
            mode=mode,
            config=self.config,
            num_samples=num_samples,
            scheduler=scheduler,
            search_alg=search_alg,
            progress_reporter=self.reporter,
            name=name,
            resources_per_trial=resources_per_trial,
            **kwargs,
        )
        best_config = analysis.best_config
        print("Best hyperparameters found were: ", best_config)
        # Get the checkpoint path of the best trial of the experiment
        model_config = {}
        trainer_config = {}
        plan_config = {}
        for key in best_config:
            if key in self.model_hyperparams:
                model_config[key] = best_config[key]
            elif key in self.trainer_hyperparams:
                trainer_config[key] = best_config[key]
            elif key in self.plan_hyperparams:
                plan_config[key] = best_config[key]
        best_checkpoint = analysis.best_checkpoint
        best_model = self.model(self.adata, **model_config)
        best_model.load(adata=self.adata, dir_path=best_checkpoint + "checkpoint")
        return best_model, analysis
