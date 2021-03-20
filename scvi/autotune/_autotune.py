import os

from ray import tune
from ray.tune import CLIReporter
from ray.tune.integration.pytorch_lightning import TuneReportCallback


class Autotune:
    """

    Hyperparameter tuning for SCVI using Ray Tune.

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
        adata,
        model,
        training_metrics: list = ["elbo_validation"],
        metric_functions: dict = {},
        model_hyperparams: dict = {},
        trainer_hyperparams: dict = {},
        plan_hyperparams: dict = {},
        num_epochs: int = 2,
    ):
        self.adata = adata
        self.model = model
        self.training_metrics = training_metrics
        self.metric_functions = metric_functions
        self.model_hyperparams = model_hyperparams
        self.trainer_hyperparams = trainer_hyperparams
        self.plan_hyperparams = plan_hyperparams
        self.metrics = training_metrics + list(self.metric_functions.keys())
        self.reporter = CLIReporter(metric_columns=self.metrics)
        self.config = {}
        for d in [model_hyperparams, trainer_hyperparams, plan_hyperparams]:
            if d is not None:
                self.config.update(d)
        self.num_epochs = num_epochs

    def _scvi_trainable(self, config, checkpoint_dir=None):
        model_config = {}
        trainer_config = {}
        plan_config = {}
        for key in config:
            if key in self.trainer_hyperparams:
                model_config[key] = config[key]
            elif key in self.trainer_hyperparams:
                trainer_config[key] = config[key]
            elif key in self.trainer_hyperparams:
                plan_config[key] = config[key]

        if checkpoint_dir:
            _model = self.model.load(
                os.path.join(checkpoint_dir, "checkpoint"), self.adata
            )
        else:
            _model = self.model(self.adata, **model_config)

        _model.train(
            **trainer_config,
            plan_kwargs=plan_config,
            callbacks=[TuneReportCallback(metrics=self.metrics, on="validation_end")],
            check_val_every_n_epoch=1,
            max_epochs=self.num_epochs,
        )
        for m in self.metric_functions:
            f = self.metric_functions[m]
            tune.report(**{m: f(_model)})
        if checkpoint_dir:
            _model.save(os.path.join(checkpoint_dir, "checkpoint"))

    def run(
        self,
        metric,
        scheduler,
        mode="min",
        name="scvi-experiment",
        num_samples=10,
        **kwargs,
    ):
        """
        Run hyper parameter tuning experiment.

        Parameters
        ----------
        metric
            Metric to optimize over in self.metrics or from self.training_funcs
        scheduler
            Ray tune scheduler for trials (e.g. ASHA).
        mode
            "min" or "max" to maximize or minimize the objective metric
        name
            Name of this experiment.
        num_samples
            Number of times to sample hyperparameters from the configuration space.

        """
        analysis = tune.run(
            self._scvi_trainable,
            metric=metric,
            mode=mode,
            config=self.config,
            num_samples=num_samples,
            scheduler=scheduler,
            progress_reporter=self.reporter,
            name=name,
        )
        best_config = analysis.best_config
        print("Best hyperparameters found were: ", best_config)
        return best_config
