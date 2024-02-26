from __future__ import annotations

from typing import Literal

from ray.tune import Tuner

from scvi._types import AnnOrMuData
from scvi.autotune._experiment import Experiment
from scvi.autotune._utils import add_rich_table_columns
from scvi.model.base import BaseModelClass


class ModelTuner:
    """``BETA`` Automated and scalable hyperparameter tuning for scvi-tools models.

    Wraps a :class:`~ray.tune.Tuner` instance attached to a scvi-tools model class.

    Parameters
    ----------
    model_cls
        A model class on which to tune hyperparameters. Must have a class property
        `_tunables` that defines tunable elements.

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
        if hasattr(self, "_model_cls"):
            raise AttributeError("Cannot reassign `model_cls`.")
        elif not hasattr(value, "_tunables"):
            raise NotImplementedError(
                f"{value} is unsupported. Please implement a `_tunables` class property to define"
                "tunable hyperparameters."
            )
        self._model_cls = value

    @property
    def registry(self) -> dict:
        from scvi.autotune._utils import get_tunables

        if hasattr(self, "_registry"):
            return self._registry

        self._registry = get_tunables(self._model_cls)
        return self._registry

    def view_registry(self, extended_info: bool = False) -> None:
        """Print the registry of tunable hyperparameters.

        Parameters
        ----------
        extended_info
            If ``True``, print extended information about the tunable hyperparameters.
        """
        import sys

        from rich.console import Console
        from rich.table import Table

        in_colab = "google.colab" in sys.modules
        force_jupyter = None if not in_colab else True
        console = Console(force_jupyter=force_jupyter)

        table = add_rich_table_columns(
            Table(title=f"{self._model_cls.__name__} tunable hyperparameters"),
            columns=["Parameter", "Default value", "Source"],
        )
        for param, metadata in self.registry.items():
            table.add_row(str(param), str(metadata["default_value"]), str(metadata["source"]))
        console.print(table)

        if extended_info:
            info_table = add_rich_table_columns(
                Table(title="Extended information"),
                columns=["Parameter", "Annotation", "Tunable type"],
            )
            for param, metadata in self.registry.items():
                info_table.add_row(
                    str(param), str(metadata["annotation"]), str(metadata["tunable_type"])
                )
            console.print(info_table)

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

    def _validate_search_space(self, search_space: dict[str, callable]) -> None:
        for param in search_space:
            if param in self.registry:
                continue
            raise ValueError(f"Parameter {param} is not in the registry.")

    def _parse_search_space(
        self,
        search_space: dict[str, callable],
    ) -> tuple[dict[str, callable], dict[str, callable]]:
        model_kwargs = {}
        train_kwargs = {}
        plan_kwargs = {}

        for param, value in search_space.items():
            type_ = self.registry[param]["tunable_type"]
            if type_ == "model":
                model_kwargs[param] = value
            elif type_ == "train":
                train_kwargs[param] = value
            elif type_ == "plan":
                plan_kwargs[param] = value

        train_kwargs["plan_kwargs"] = plan_kwargs
        return model_kwargs, train_kwargs

    def _get_trainable(self, experiment: Experiment) -> callable:
        import os

        from lightning.pytorch.callbacks import Callback
        from lightning.pytorch.loggers import TensorBoardLogger
        from ray.train import get_context
        from ray.tune import with_parameters, with_resources
        from ray.tune.integration.pytorch_lightning import TuneReportCheckpointCallback

        from scvi import settings

        tune_callback_cls = type(
            "_TuneReportCallback", (TuneReportCheckpointCallback, Callback), {}
        )

        def _trainable(search_space: dict[str, callable], experiment: Experiment) -> None:
            settings.seed = experiment.seed

            callbacks = [
                tune_callback_cls(
                    metrics=experiment.metrics,
                    on="validation_end",
                    save_checkpoints=False,
                )
            ]

            base_dir = os.path.join(experiment.logging_dir, experiment.experiment_name)
            log_dir = os.path.join(base_dir, f"{get_context().get_trial_name()}_tensorboard")
            logger = TensorBoardLogger(log_dir)

            model_kwargs, train_kwargs = self._parse_search_space(search_space)
            model_kwargs = {**experiment.model_kwargs, **model_kwargs}
            train_kwargs = {
                "max_epochs": experiment.max_epochs,
                "accelerator": experiment.accelerator,
                "devices": experiment.devices,
                "check_val_every_n_epoch": 1,
                "enable_progress_bar": False,
                "logger": logger,
                "callbacks": callbacks,
                **experiment.train_kwargs,
                **train_kwargs,
            }

            getattr(self.model_cls, experiment.setup_method_name)(
                experiment.adata,
                **experiment.setup_args,
            )
            model = self.model_cls(experiment.adata, **model_kwargs)
            model.train(**train_kwargs)

        trainable_with_params = with_parameters(_trainable, experiment=experiment)
        return with_resources(trainable_with_params, resources=experiment.resources)

    def _configure_experiment(
        self,
        adata: AnnOrMuData,
        metrics: str | list[str],
        mode: Literal["min", "max"],
        search_space: dict[str, callable],
        num_samples: int,
        model_kwargs: dict | None = None,
        train_kwargs: dict | None = None,
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
        self._validate_search_space(search_space)
        experiment = Experiment(
            self.model_cls,
            adata,
            metrics,
            mode,
            search_space,
            num_samples,
            model_kwargs=model_kwargs,
            train_kwargs=train_kwargs,
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
        model_kwargs: dict | None = None,
        train_kwargs: dict | None = None,
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
            search specifications.
        num_samples
            Total number of hyperparameter configurations to sample from the search space. Passed
            into :class:`~ray.tune.tune_config.TuneConfig`.
        model_kwargs
            Keyword arguments passed to the model class's constructor. Arguments must not overlap
            with those in ``search_space``.
        train_kwargs
            Keyword arguments passed to the model's ``train`` method. Arguments must not overlap
            with those in ``search_space``.
        max_epochs
            Maximum number of epochs to train each hyperparameter configuration.
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
            Random seed to use for the experiment.
        resources
            Dictionary of resources to allocate for each trial in the experiment. Available keys
            include:

            * ``"cpu"``: number of CPU cores
            * ``"gpu"``: number of GPUs
            * ``"memory"``: amount of memory

            If not provided, defaults to using all available resources. Note that fractional
            allocations are supported.
        experiment_name
            Name of the experiment, used for logging purposes. Defaults to a unique
            string formatted with the current timestamp and model class name.
        logging_dir
            Directory to store experiment logs. Defaults to a directory named ``"autotune"``
            in :attr:``scvi.settings.logging_dir``.
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
            model_kwargs=model_kwargs,
            train_kwargs=train_kwargs,
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
        result_grid = self._configure_tuner(experiment).fit()
        experiment.result_grid = result_grid
        self._update_history(experiment)

        return experiment

    def __repr__(self) -> str:
        return f"ModelTuner for {self.manager.model_cls.__name__}"
