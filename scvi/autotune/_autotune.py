from ray import tune
from ray.tune import CLIReporter
from ray.tune.integration.pytorch_lightning import TuneReportCallback
from ray.tune.schedulers import ASHAScheduler


class Autotune:
    def __init__(self, adata, model, name, metric, metrics, config, num_samples=10):
        self.adata = adata
        self.model = model
        self.name = name
        self.metric = metric
        self.metrics = metrics
        self.config = config
        self.num_samples = num_samples
        self.scheduler = ASHAScheduler(max_t=10, grace_period=1, reduction_factor=2)
        self.reporter = CLIReporter(metric_columns=self.metrics)

    def set_reporter(self, reporter):
        self.reporter = reporter

    def set_scheduler(self, scheduler):
        self.scheduler = scheduler

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

    def __call__(self):
        analysis = tune.run(
            tune.with_parameters(
                self.train_model_tune, adata=self.adata, metrics=self.metrics
            ),
            metric=self.metric,
            mode="min",
            config=self.config,
            num_samples=self.num_samples,
            scheduler=self.scheduler,
            progress_reporter=self.reporter,
            name=self.name,
        )
        best_config = analysis.best_config
        print("Best hyperparameters found were: ", best_config)

        return self.model(self.adata, **best_config)
