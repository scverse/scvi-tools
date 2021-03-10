from ray import tune
from ray.tune import CLIReporter
from ray.tune.integration.pytorch_lightning import TuneReportCallback
from ray.tune.schedulers import ASHAScheduler


class Autotune:
    def __init__(self, adata, model, n_samples=10):
        self.adata = adata
        self.model = model
        self.n_samples = n_samples

    def train_scvi_tune(self, config, adata=None):
        _model = self.model(adata, **config)
        callback = TuneReportCallback(
            {"loss": "avg_val_loss", "mean_accuracy": "avg_val_accuracy"},
            on="validation_end",
        )
        _model.train(callbacks=[callback], max_epochs=2)

    def tune_scvi_asha(self):
        config = {
            "dropout_rate": tune.loguniform(1e-4, 1e-1),
        }
        scheduler = ASHAScheduler(max_t=10, grace_period=1, reduction_factor=2)
        reporter = CLIReporter(
            metric_columns=["loss", "mean_accuracy", "training_iteration"]
        )

        analysis = tune.run(
            tune.with_parameters(self.train_scvi_tune, adata=self.adata),
            metric="loss",
            mode="min",
            config=config,
            num_samples=self.n_samples,
            scheduler=scheduler,
            progress_reporter=reporter,
            name="tune_scvi_asha",
        )

        print("Best hyperparameters found were: ", analysis.best_config)
