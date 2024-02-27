from __future__ import annotations

from typing import Literal

from ray.tune import Tuner

from scvi._types import AnnOrMuData
from scvi.autotune._experiment import Experiment
from scvi.model.base import BaseModelClass


class ModelTuner:
    """``BETA`` Automated and scalable hyperparameter tuning for scvi-tools models.

    Wraps a :class:`~ray.tune.Tuner` instance attached to a scvi-tools model class.

    Parameters
    ----------
    model_cls
        A model class on which to tune hyperparameters.

    Examples
    --------
    >>> import anndata
    >>> import scvi
    >>> adata = anndata.read_h5ad(path_to_h5ad)
    >>> model_cls = scvi.model.SCVI
    >>> model_cls.setup_anndata(adata)
    >>> tuner = scvi.autotune.ModelTuner(model_cls)
    >>> results = tuner.fit(adata, metric="validation_loss")

    Notes
    -----
    See further usage examples in the following tutorials:

    1. :doc:`/tutorials/notebooks/tuning/autotune_new_model`
    2. :doc:`/tutorials/notebooks/tuning/autotune_scvi`

    Lifecycle: beta
    """

    def __init__(self, model_cls: BaseModelClass):
        self.model_cls = model_cls

    @property
    def model_cls(self) -> BaseModelClass:
        """The model class associated with the manager."""
        return self._model_cls

    @model_cls.setter
    def model_cls(self, value: BaseModelClass) -> None:
        """Assigns the model class.

        Checks if the model class is supported.

        Raises
        ------
        AttributeError
            If ``model_cls`` is already assigned.
        NotImplementedError
            If ``value`` is unsupported.
        """
        self._model_cls = value

    @property
    def history(self) -> dict[str, Experiment]:
        """History of tuning experiments.

        Can be cleared with :meth:`~scvi.autotune.TunerManager.clear_history`.
        """
        if hasattr(self, "_history"):
            return self._history

        self._history = {}
        return self._history

    def _update_history(self, experiment: Experiment) -> None:
        """Add an experiment to the history.

        Parameters
        ----------
        experiment
            The experiment to add to the history.
        """
        history = self.history
        history[experiment.id] = experiment
        self._history = history

    def clear_history(self, key: str | None = None) -> None:
        """Clear the history of tuning experiments.

        Parameters
        ----------
        key
            If provided, only the experiment with the given key will be removed. Otherwise, the
            entire history will be cleared.

        Raises
        ------
        KeyError
            If ``key`` is not found in the history.
        """
        if key is None:
            self._history = {}
        elif key not in self.history:
            raise KeyError(f"Experiment with key {key} not found in history.")
        else:
            history = self.history
            del history[key]
            self._history = history

    def _get_trainable(self, experiment: Experiment) -> callable:
        from os.path import join

        from lightning.pytorch.callbacks import Callback
        from lightning.pytorch.loggers import TensorBoardLogger
        from ray.train import get_context
        from ray.tune import with_parameters, with_resources
        from ray.tune.integration.pytorch_lightning import TuneReportCheckpointCallback

        tune_callback_cls = type(
            "_TuneReportCallback", (TuneReportCheckpointCallback, Callback), {}
        )

        def _trainable(search_space: dict[str, callable], experiment: Experiment) -> None:
            from scvi import settings

            settings.seed = experiment.seed

            callbacks = [
                tune_callback_cls(
                    metrics=experiment.metrics,
                    on="validation_end",
                    save_checkpoints=False,
                )
            ]
            log_dir = join(
                experiment.logging_dir,
                experiment.experiment_name,
                f"{get_context().get_trial_name()}_tensorboard",
            )

            model_args = search_space.get("model_kwargs", {})
            train_args = search_space.get("train_kwargs", {})
            train_args = {
                "max_epochs": experiment.max_epochs,
                "accelerator": experiment.accelerator,
                "devices": experiment.devices,
                "check_val_every_n_epoch": 1,
                "enable_progress_bar": False,
                "logger": TensorBoardLogger(log_dir),
                "callbacks": callbacks,
                **train_args,
            }

            getattr(self.model_cls, experiment.setup_method_name)(
                experiment.adata,
                **experiment.setup_args,
            )
            model = self.model_cls(experiment.adata, **model_args)
            model.train(**train_args)

        trainable_with_params = with_parameters(_trainable, experiment=experiment)
        return with_resources(trainable_with_params, resources=experiment.resources)

    def _configure_experiment(
        self,
        adata: AnnOrMuData,
        metrics: str | list[str],
        mode: Literal["min", "max"],
        search_space: dict[str, callable],
        num_samples: int,
        max_epochs: int | None = None,
        scheduler: str | None = None,
        searcher: str | None = None,
        seed: int | None = None,
        resources: dict[str, float] | None = None,
        experiment_name: str | None = None,
        logging_dir: str | None = None,
        scheduler_kwargs: dict | None = None,
        searcher_kwargs: dict | None = None,
    ) -> Experiment:
        experiment = Experiment(
            self.model_cls,
            adata,
            metrics,
            mode,
            search_space,
            num_samples,
            scheduler=scheduler,
            searcher=searcher,
            seed=seed,
            resources=resources,
            experiment_name=experiment_name,
            logging_dir=logging_dir,
            scheduler_kwargs=scheduler_kwargs,
            searcher_kwargs=searcher_kwargs,
        )
        return experiment

    def _configure_tuner(self, experiment: Experiment) -> Tuner:
        from ray.train import RunConfig
        from ray.tune.tune_config import TuneConfig

        trainable = self._get_trainable(experiment)
        tune_config = TuneConfig(
            scheduler=experiment.scheduler,
            search_alg=experiment.searcher,
            num_samples=experiment.num_samples,
        )
        run_config = RunConfig(
            name=experiment.experiment_name,
            storage_path=experiment.logging_dir,
            local_dir=experiment.logging_dir,
            log_to_file=True,
            verbose=1,
        )
        return Tuner(
            trainable=trainable,
            param_space=experiment.search_space,
            tune_config=tune_config,
            run_config=run_config,
        )

    def fit(
        self,
        adata: AnnOrMuData,
        metrics: str | list[str],
        mode: Literal["min", "max"],
        search_space: dict[str, callable],
        num_samples: int,
        max_epochs: int | None = None,
        scheduler: Literal["asha", "hyperband", "median", "pbt", "fifo"] = "asha",
        searcher: Literal["hyperopt", "random"] = "hyperopt",
        seed: int | None = None,
        resources: dict[str, float] | None = None,
        experiment_name: str | None = None,
        logging_dir: str | None = None,
        scheduler_kwargs: dict | None = None,
        searcher_kwargs: dict | None = None,
    ) -> Experiment:
        """Run a specified hyperparameter sweep for the associated model class.

        Parameters
        ----------
        adata
            :class:`~anndata.AnnData` or :class:`~mudata.MuData` that has been setup with the
            associated model class.
        metrics
            Either a single metric or a list of metrics to track during the experiment. If a list
            is provided, the primary metric will be the first element in the list.
        mode
            Optimization mode for the primary metric. One of ``"min"`` or ``"max"``.
        search_space
            Dictionary of hyperparameter names and their respective search spaces. See
            the `API <https://docs.ray.io/en/latest/tune/api/search_space.html>`_ for available
            search specifications. Must only contain the following top-level keys:

            * ``"model_args"``: parameters to pass to the model constructor.
            * ``"train_args"``: parameters to pass to the model's ``train`` method.
        num_samples
            Total number of hyperparameter configurations to sample from the search space. Passed
            into :class:`~ray.tune.tune_config.TuneConfig`.
        max_epochs
            Maximum number of epochs to train hyperparameter configurations.
        scheduler
            Ray Tune scheduler to use. One of the following:

            * ``"asha"``: :class:`~ray.tune.schedulers.AsyncHyperBandScheduler`
            * ``"hyperband"``: :class:`~ray.tune.schedulers.HyperBandScheduler`
            * ``"median"``: :class:`~ray.tune.schedulers.MedianStoppingRule`
            * ``"pbt"``: :class:`~ray.tune.schedulers.PopulationBasedTraining`
            * ``"fifo"``: :class:`~ray.tune.schedulers.FIFOScheduler`

            Note that that not all schedulers are compatible with all search algorithms. See the
            `documentation <https://docs.ray.io/en/latest/tune/key-concepts.html#schedulers>`_
            for more details.
        searcher
            Ray Tune search algorithm to use. One of the following:

            * ``"hyperopt"``: :class:`~ray.tune.search.hyperopt.HyperOptSearch`
            * ``"random"``: :class:`~ray.tune.search.basic_variant.BasicVariantGenerator`
        seed
            Random seed to use for the experiment. Propagated to :attr:`~scvi.settings.seed` and
            search algorithms. If not provided, defaults to :attr:`~scvi.settings.seed`.
        resources
            Dictionary of resources to allocate per trial in the experiment. Available keys
            include:

            * ``"cpu"``: number of CPU cores
            * ``"gpu"``: number of GPUs
            * ``"memory"``: amount of memory

            If not provided, defaults to using all available resources. Note that fractional
            allocations are supported.
        experiment_name
            Name of the experiment, used for logging purposes. Defaults to a unique ID.
        logging_dir
            Base directory to store experiment logs. Defaults to :attr:``scvi.settings.logging_dir``.
        scheduler_kwargs
            Additional keyword arguments to pass to the scheduler.
        searcher_kwargs
            Additional keyword arguments to pass to the search algorithm.
        """
        experiment = self._configure_experiment(
            adata,
            metrics,
            mode,
            search_space,
            num_samples,
            max_epochs=max_epochs,
            scheduler=scheduler,
            searcher=searcher,
            seed=seed,
            resources=resources,
            experiment_name=experiment_name,
            logging_dir=logging_dir,
            scheduler_kwargs=scheduler_kwargs,
            searcher_kwargs=searcher_kwargs,
        )
        experiment.result_grid = self._configure_tuner(experiment).fit()
        self._update_history(experiment)
        return experiment

    def __repr__(self) -> str:
        return f"ModelTuner for {self.manager.model_cls.__name__}"
