from __future__ import annotations

import logging
from typing import Literal

from ray.tune import Tuner

from scvi._types import AnnOrMuData
from scvi.autotune._experiment import Experiment
from scvi.autotune._utils import add_rich_table_columns
from scvi.model.base import BaseModelClass

logger = logging.getLogger(__name__)


class TunerManager:
    """Internal manager for validation of inputs from :class:`~scvi.autotune.ModelTuner`.

    Keeps a history of the tuning experiments and provides a view of the tunable hyperparameters.

    Parameters
    ----------
    model_cls
        A model class on which to tune hyperparameters. Must implement a ``_tunables`` class
        property that defines tunable parameters.
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
                    metrics=[experiment.metric] + experiment.additional_metrics,
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

    def configure_tuner(
        self,
        adata: AnnOrMuData,
        metric: str,
        mode: Literal["min", "max"],
        search_space: dict[str, callable],
        num_samples: int,
        additional_metrics: list[str] | None = None,
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
    ) -> tuple[Tuner, Experiment]:
        import os

        from ray.air.config import RunConfig
        from ray.tune.tune_config import TuneConfig

        self._validate_search_space(search_space)
        experiment = Experiment(
            self.model_cls,
            adata,
            metric,
            mode,
            search_space,
            num_samples,
            additional_metrics=additional_metrics,
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
        self._update_history(experiment)

        trainable = self._get_trainable(experiment)
        tune_config = TuneConfig(
            scheduler=experiment.scheduler,
            search_alg=experiment.searcher,
            num_samples=experiment.num_samples,
        )
        save_dir = os.path.join(experiment.logging_dir)
        run_config = RunConfig(
            name=experiment.experiment_name,
            storage_path=save_dir,
            local_dir=save_dir,
            log_to_file=True,
            verbose=1,
        )
        tuner = Tuner(
            trainable=trainable,
            param_space=experiment.search_space,
            tune_config=tune_config,
            run_config=run_config,
        )
        return tuner, experiment
