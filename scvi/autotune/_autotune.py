from typing import Optional

from ray import tune
from ray.tune import CLIReporter
from ray.tune.integration.pytorch_lightning import TuneReportCallback
from ray.tune.schedulers import ASHAScheduler

"""
user needs to be able to set search space for model kwargs, trainer kwargs, training plan kwargs. 
track validation metrics, as well as metrics after training that come from metric functions
"""


class Autotune:
    def __init__(
        self,
        adata,
        model,
        name: str = "scvi-experiment",
        metric: str = "",
        training_metrics: dict = None,
        post_training_metric_functions: dict = None,
        model_hyperparams: dict = None,
        trainer_hyperparams: dict = None,
        training_plan_hyperparams: dict = None,
        search_alg: Searcher = None,
        scheduler: TrialScheduler = None,
        checkpoint_dir: str = "experiment-checkpoints"
    ):
        if not training_metrics
        
        
    def train_model_tune(self, config, adata, metrics):
        _current_model = self.model(adata, **config)
        _current_model.train(
            check_val_every_n_epoch=1,
            callbacks=[
                TuneReportCallback(
                    metrics,
                    on="validation_end",
                )
            ],
            max_epochs=2,  # remove this later
        )

    def run(self, metric, mode, name, num_samples):
        analysis = tune.run(
            tune.with_parameters(
                self.train_model_tune, adata=self.adata, metrics=self.metrics
            ),
            metric=metric,
            mode="min",
            config=self.config,
            num_samples=num_samples,
            scheduler=self.scheduler,
            progress_reporter=self.reporter,
            name=name,
        )
        best_config = analysis.best_config
        print("Best hyperparameters found were: ", best_config)

        return self.model(self.adata, **best_config)
