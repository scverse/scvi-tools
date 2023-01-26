import inspect
import logging
import os
import warnings
from collections import OrderedDict
from datetime import datetime
from typing import Any, Callable, List, Optional, Tuple

import rich
from chex import dataclass

try:
    import ray
    from ray import air, tune
    from ray.tune.integration.pytorch_lightning import TuneReportCallback
except ImportError:
    pass

from scvi._decorators import dependencies
from scvi._types import AnnOrMuData
from scvi.data._constants import _SETUP_ARGS_KEY, _SETUP_METHOD_NAME
from scvi.model.base import BaseModelClass
from scvi.utils import InvalidParameterError

from ._defaults import COLORS, COLUMN_KWARGS, DEFAULTS, TUNABLE_TYPES
from ._types import TunableMeta
from ._utils import in_notebook

logger = logging.getLogger(__name__)


@dataclass
class TuneAnalysis:
    """Dataclass for storing results from a tuning experiment."""

    model_kwargs: dict
    train_kwargs: dict
    metric: float
    additional_metrics: dict
    search_space: dict
    results: Any


class TunerManager:
    """
    Internal manager for validation of inputs from :class:`~scvi.autotune.ModelTuner`.

    Parameters
    ----------
    model_cls
        A model class on which to tune hyperparameters. Must have a class property
        `_tunables` that defines tunable elements.
    """

    def __init__(self, model_cls: BaseModelClass):
        self._model_cls = self._validate_model_cls(model_cls)
        self._defaults = self._get_defaults(self._model_cls)
        self._registry = self._get_registry(self._model_cls)

    @staticmethod
    def _validate_model_cls(model_cls: BaseModelClass) -> BaseModelClass:
        """Checks if the model class is supported."""
        if not hasattr(model_cls, "_tunables"):
            raise NotImplementedError(
                f"{model_cls} is unsupported. Please implement a `_tunables` class "
                "property to define tunable hyperparameters."
            )
        return model_cls

    @staticmethod
    def _get_defaults(model_cls: BaseModelClass) -> dict:
        """Returns the model class's default search space if available."""
        if model_cls not in DEFAULTS:
            warnings.warn(
                f"No default search space available for {model_cls.__name__}.",
                UserWarning,
            )
        return DEFAULTS.get(model_cls, {})

    @staticmethod
    def _get_registry(model_cls: BaseModelClass) -> dict:
        """
        Returns the model class's registry of tunable hyperparameters and metrics.

        For a given model class, checks whether a `_tunables` class property has been
        defined. If so, iterates through the attribute to recursively discover tunable
        hyperparameters.

        Parameters
        ----------
        model_cls
            A validated :class:`~scvi.model.base.BaseModelClass`.

        Returns
        -------
        registry: dict
            A dictionary with the following keys:

            * ``"tunables"``: a dictionary of tunable hyperparameters and metadata
            * ``"metrics"``: a dictionary of available metrics and metadata
        """

        def _cls_to_tunable_type(cls: Any) -> str:
            for tunable_type, cls_list in TUNABLE_TYPES.items():
                if any([issubclass(cls, c) for c in cls_list]):
                    return tunable_type
            return None

        def _parse_func_params(func: Callable, parent: Any, tunable_type: str) -> dict:
            # get function kwargs that are tunable
            tunables = {}
            for param, metadata in inspect.signature(func).parameters.items():
                if not isinstance(metadata.annotation, TunableMeta):
                    continue

                default_val = None
                if metadata.default is not inspect.Parameter.empty:
                    default_val = metadata.default

                annotation = metadata.annotation.__args__[0]
                if hasattr(annotation, "__args__"):
                    # e.g. if type is Literal, get its arguments
                    annotation = annotation.__args__
                else:
                    annotation = annotation.__name__

                tunables[param] = {
                    "tunable_type": tunable_type,
                    "default_value": default_val,
                    "source": parent.__name__,
                    "annotation": annotation,
                }
            return tunables

        def _get_tunables(
            attr: Any, parent: Any = None, tunable_type: Optional[str] = None
        ) -> dict:
            tunables = {}
            if inspect.isfunction(attr):
                return _parse_func_params(attr, parent, tunable_type)
            for child in getattr(attr, "_tunables", []):
                tunables.update(
                    _get_tunables(
                        child, parent=attr, tunable_type=_cls_to_tunable_type(attr)
                    )
                )
            return tunables

        def _get_metrics(model_cls: BaseModelClass) -> OrderedDict:
            # TODO: discover more metrics
            return {"validation_loss": "min"}

        registry = {
            "tunables": _get_tunables(model_cls),
            "metrics": _get_metrics(model_cls),
        }
        return registry

    def _get_search_space(self, search_space: dict) -> Tuple[dict, dict]:
        """Parses a compact search space into separate kwargs dictionaries."""
        model_kwargs = {}
        train_kwargs = {}
        plan_kwargs = {}
        tunables = self._registry["tunables"]

        for param, value in search_space.items():
            _type = tunables[param]["tunable_type"]
            if _type == "model":
                model_kwargs[param] = value
            elif _type == "train":
                train_kwargs[param] = value
            elif _type == "plan":
                plan_kwargs[param] = value

        train_kwargs["plan_kwargs"] = plan_kwargs
        return model_kwargs, train_kwargs

    @dependencies("ray.tune")
    def _validate_search_space(self, search_space: dict, use_defaults: bool) -> dict:
        """Validates a search space against the hyperparameter registry."""
        # validate user-provided search space
        for param in search_space:
            if param in self._registry["tunables"]:
                continue
            raise ValueError(
                f"Provided parameter `{param}` is invalid for {self._model_cls.__name__}."
                " Please see available parameters with `ModelTuner.info()`. "
            )

        # add defaults if requested
        _search_space = {}
        if use_defaults:

            # parse defaults into tune sample functions
            for param, metadata in self._defaults.items():
                sample_fn = getattr(tune, metadata["fn"])
                fn_args = metadata.get("args", [])
                fn_kwargs = metadata.get("kwargs", {})
                _search_space[param] = sample_fn(*fn_args, **fn_kwargs)

            # exclude defaults if requested
            logger.info(
                f"Merging search space with defaults for {self._model_cls.__name__}."
            )

        # priority given to user-provided search space
        _search_space.update(search_space)
        return _search_space

    def _validate_metrics(
        self, metric: str, additional_metrics: List[str]
    ) -> OrderedDict:
        """Validates a set of metrics against the metric registry."""
        registry_metrics = self._registry["metrics"]
        _metrics = OrderedDict()

        # validate primary metric
        if metric not in registry_metrics:
            raise ValueError(
                f"Provided metric `{metric}` is invalid for {self._model_cls.__name__}. "
                "Please see available metrics with `ModelTuner.info()`. ",
            )
        _metrics[metric] = registry_metrics[metric]

        # validate additional metrics
        for m in additional_metrics:
            if m not in registry_metrics:
                warnings.warn(
                    f"Provided metric {m} is invalid for {self._model_cls.__name__}. "
                    "Please see available metrics with `ModelTuner.info()`. "
                    "Ignoring metric.",
                    UserWarning,
                )
                continue
            _metrics[m] = registry_metrics[m]

        return _metrics

    @staticmethod
    def _get_primary_metric_and_mode(metrics: OrderedDict) -> Tuple[str, str]:
        metric = list(metrics.keys())[0]
        mode = metrics[metric]
        return metric, mode

    @dependencies("ray.tune")
    def _validate_scheduler(
        self, scheduler: str, metrics: OrderedDict, scheduler_kwargs: dict
    ) -> Any:
        """Validates a trial scheduler."""
        metric, mode = self._get_primary_metric_and_mode(metrics)
        _kwargs = {"metric": metric, "mode": mode}

        if scheduler == "asha":
            _default_kwargs = {
                "max_t": 100,
                "grace_period": 1,
                "reduction_factor": 2,
            }
            _scheduler = tune.schedulers.AsyncHyperBandScheduler
        elif scheduler == "hyperband":
            _default_kwargs = {}
            _scheduler = tune.schedulers.HyperBandScheduler
        elif scheduler == "median":
            _default_kwargs = {}
            _scheduler = tune.schedulers.MedianStoppingRule
        elif scheduler == "pbt":
            _default_kwargs = {}
            _scheduler = tune.schedulers.PopulationBasedTraining
        elif scheduler == "fifo":
            _default_kwargs = {}
            _scheduler = tune.schedulers.FIFOScheduler

        # prority given to user-provided scheduler kwargs
        _default_kwargs.update(scheduler_kwargs)
        _kwargs.update(_default_kwargs)
        return _scheduler(**_kwargs)

    @dependencies(["ray.tune", "hyperopt"])
    def _validate_search_algorithm(
        self, searcher: str, metrics: OrderedDict, searcher_kwargs: dict
    ) -> Any:
        """Validates a hyperparameter search algorithm."""
        metric, mode = self._get_primary_metric_and_mode(metrics)

        if searcher == "random":
            _default_kwargs = {}
            _searcher = tune.search.basic_variant.BasicVariantGenerator
        elif searcher == "hyperopt":
            _default_kwargs = {
                "metric": metric,
                "mode": mode,
            }
            # tune does not import hyperopt by default
            tune.search.SEARCH_ALG_IMPORT[searcher]()
            _searcher = tune.search.hyperopt.HyperOptSearch

        # prority given to user-provided searcher kwargs
        _default_kwargs.update(searcher_kwargs)
        return _searcher(**_default_kwargs)

    def _validate_scheduler_and_search_algorithm(
        self,
        scheduler: str,
        searcher: str,
        metrics: OrderedDict,
        scheduler_kwargs: dict,
        searcher_kwargs: dict,
    ) -> Tuple[Any, Any]:
        """Validates a scheduler and search algorithm pair for compatibility."""
        supported = ["asha", "hyperband", "median", "pbt", "fifo"]
        if scheduler not in supported:
            raise InvalidParameterError("scheduler", scheduler, supported)

        supported = ["random", "hyperopt"]
        if searcher not in supported:
            raise InvalidParameterError("searcher", searcher, supported)

        if scheduler not in ["asha", "median", "hyperband"] and searcher not in [
            "random",
            "grid",
        ]:
            raise ValueError(
                f"Provided scheduler {scheduler} is incompatible with the provided "
                f"searcher {searcher}. Please see "
                "https://docs.ray.io/en/latest/tune/key-concepts.html#schedulers for more info."
            )

        _scheduler = self._validate_scheduler(scheduler, metrics, scheduler_kwargs)
        _searcher = self._validate_search_algorithm(searcher, metrics, searcher_kwargs)
        return _scheduler, _searcher

    @staticmethod
    @dependencies("ray.tune")
    def _validate_reporter(
        reporter: bool, search_space: dict, metrics: OrderedDict
    ) -> Any:
        """Validates a reporter depending on the execution environment."""
        _metric_keys = list(metrics.keys())
        _param_keys = list(search_space.keys())
        _kwargs = {
            "metric_columns": _metric_keys,
            "parameter_columns": _param_keys,
            "metric": _metric_keys[0],
            "mode": metrics[_metric_keys[0]],
        }

        if not reporter:
            _reporter = None
        elif in_notebook():
            _reporter = tune.JupyterNotebookReporter(**_kwargs)
        else:
            _reporter = tune.CLIReporter(**_kwargs)

        return _reporter

    def _validate_resources(self, resources: dict) -> dict:
        """Validates a resource-use specification."""
        # TODO: perform resource checking
        return resources

    def _get_setup_info(self, adata: AnnOrMuData) -> Tuple[str, dict]:
        """Retrieves the method and kwargs used for setting up `adata` with the model class."""
        manager = self._model_cls._get_most_recent_anndata_manager(adata)
        setup_method_name = manager._registry.get(_SETUP_METHOD_NAME, "setup_anndata")
        setup_args = manager._get_setup_method_args().get(_SETUP_ARGS_KEY, {})
        return setup_method_name, setup_args

    @dependencies("ray.tune")
    def _get_trainable(
        self,
        adata: AnnOrMuData,
        metrics: OrderedDict,
        resources: dict,
        setup_method_name: str,
        setup_kwargs: dict,
        max_epochs: int,
    ) -> Callable:
        """Returns a trainable function consumable by :class:`~ray.tune.Tuner`."""

        def _trainable(
            search_space: dict,
            *,
            model_cls: BaseModelClass,
            adata: AnnOrMuData,
            metric: str,
            setup_method_name: str,
            setup_kwargs: dict,
            max_epochs: int,
            use_gpu: bool,
        ) -> None:
            model_kwargs, train_kwargs = self._get_search_space(search_space)
            getattr(model_cls, setup_method_name)(adata, **setup_kwargs)
            model = model_cls(adata, **model_kwargs)
            monitor = TuneReportCallback(metric, on="validation_end")
            model.train(
                max_epochs=max_epochs,
                use_gpu=use_gpu,
                check_val_every_n_epoch=1,
                callbacks=[monitor],
                enable_progress_bar=False,
                **train_kwargs,
            )

        _wrap_params = tune.with_parameters(
            _trainable,
            model_cls=self._model_cls,
            adata=adata,
            metric=list(metrics.keys())[0],
            setup_method_name=setup_method_name,
            setup_kwargs=setup_kwargs,
            max_epochs=max_epochs,
            use_gpu=resources.get("gpu", 0) > 0,
        )
        return tune.with_resources(_wrap_params, resources=resources)

    def _validate_experiment_name_and_logging_dir(
        self, experiment_name: Optional[str], logging_dir: Optional[str]
    ) -> Tuple[str, str]:
        if experiment_name is None:
            experiment_name = "tune_"
            experiment_name += self._model_cls.__name__.lower() + "_"
            experiment_name += datetime.now().strftime("%Y-%m-%d-%H:%M:%S")
        if logging_dir is None:
            logging_dir = os.path.join(os.getcwd(), "ray")
        return experiment_name, logging_dir

    @dependencies(["ray.tune", "ray.air"])
    def _get_tuner(
        self,
        adata: AnnOrMuData,
        *,
        metric: Optional[str] = None,
        additional_metrics: Optional[List[str]] = None,
        search_space: Optional[dict] = None,
        use_defaults: bool = False,
        num_samples: Optional[int] = None,
        max_epochs: Optional[int] = None,
        scheduler: Optional[str] = None,
        scheduler_kwargs: Optional[dict] = None,
        searcher: Optional[str] = None,
        searcher_kwargs: Optional[dict] = None,
        reporter: bool = True,
        resources: Optional[dict] = None,
        experiment_name: Optional[str] = None,
        logging_dir: Optional[str] = None,
    ) -> Tuple[Any, dict]:
        metric = (
            metric or self._get_primary_metric_and_mode(self._registry["metrics"])[0]
        )
        additional_metrics = additional_metrics or []
        search_space = search_space or {}
        num_samples = num_samples or 10  # TODO: better default
        max_epochs = max_epochs or 100  # TODO: better default
        scheduler = scheduler or "asha"
        scheduler_kwargs = scheduler_kwargs or {}
        searcher = searcher or "random"
        searcher_kwargs = searcher_kwargs or {}
        resources = resources or {}

        _metrics = self._validate_metrics(metric, additional_metrics)
        _search_space = self._validate_search_space(search_space, use_defaults)
        _scheduler, _searcher = self._validate_scheduler_and_search_algorithm(
            scheduler, searcher, _metrics, scheduler_kwargs, searcher_kwargs
        )
        _reporter = self._validate_reporter(reporter, _search_space, _metrics)
        _resources = self._validate_resources(resources)
        _setup_method_name, _setup_args = self._get_setup_info(adata)
        _trainable = self._get_trainable(
            adata,
            _metrics,
            _resources,
            _setup_method_name,
            _setup_args,
            max_epochs,
        )
        _experiment_name, _logging_dir = self._validate_experiment_name_and_logging_dir(
            experiment_name, logging_dir
        )

        tune_config = tune.tune_config.TuneConfig(
            scheduler=_scheduler,
            search_alg=_searcher,
            num_samples=num_samples,
        )
        run_config = air.config.RunConfig(
            name=_experiment_name,
            local_dir=_logging_dir,
            progress_reporter=_reporter,
            log_to_file=True,
            verbose=1,
        )
        tuner = tune.Tuner(
            trainable=_trainable,
            param_space=_search_space,
            tune_config=tune_config,
            run_config=run_config,
        )
        config = {
            "metrics": _metrics,
            "search_space": _search_space,
        }
        return tuner, config

    def _get_analysis(self, results: Any, config: dict) -> TuneAnalysis:
        metrics = config["metrics"]
        search_space = config["search_space"]
        metric, mode = self._get_primary_metric_and_mode(metrics)

        result = results.get_best_result(metric=metric, mode=mode)
        model_kwargs, train_kwargs = self._get_search_space(result.config)
        metric_values = {}
        for m in metrics:
            if m == metric:
                continue
            metric_values[m] = result.metrics[m]

        return TuneAnalysis(
            model_kwargs=model_kwargs,
            train_kwargs=train_kwargs,
            metric={"metric": metric, "mode": mode, "value": result.metrics[metric]},
            additional_metrics=metric_values,
            search_space=search_space,
            results=results,
        )

    @staticmethod
    def _add_columns(table: rich.table.Table, columns: List[str]) -> rich.table.Table:
        """Adds columns to a :class:`~rich.table.Table` with default formatting."""
        for i, column in enumerate(columns):
            table.add_column(column, style=COLORS[i], **COLUMN_KWARGS)
        return table

    @staticmethod
    @dependencies("ray")
    def _get_resources(available: bool = False) -> dict:
        # TODO: need a cleaner way to do this as it starts a ray instance
        ray.init(logging_level=logging.ERROR)
        if available:
            resources = ray.available_resources()
        else:
            resources = ray.cluster_resources()
        ray.shutdown()
        return resources

    def _view_registry(
        self, show_additional_info: bool = False, show_resources: bool = False
    ) -> None:
        """
        Displays a summary of the model class's registry and available resources.

        Parameters
        ----------
        show_additional_info
            Whether to show additional information about the model class's registry.
        show_resources
            Whether to show available resources.
        """
        console = rich.console.Console(force_jupyter=in_notebook())
        console.print(f"ModelTuner registry for {self._model_cls.__name__}")

        tunables_table = self._add_columns(
            rich.table.Table(title="Tunable hyperparameters"),
            ["Hyperparameter", "Default value", "Source"],
        )
        for param, metadata in self._registry["tunables"].items():
            tunables_table.add_row(
                str(param),
                str(metadata["default_value"]),
                str(metadata["source"]),
            )
        console.print(tunables_table)

        if show_additional_info:
            additional_info_table = self._add_columns(
                rich.table.Table(title="Additional information", width=100),
                ["Hyperparameter", "Annotation", "Tunable type"],
            )
            for param, metadata in self._registry["tunables"].items():
                additional_info_table.add_row(
                    str(param),
                    str(metadata["annotation"]),
                    str(metadata["tunable_type"]),
                )
            console.print(additional_info_table)

        metrics_table = self._add_columns(
            rich.table.Table(title="Available metrics"),
            ["Metric", "Mode"],
        )
        for metric, mode in self._registry["metrics"].items():
            metrics_table.add_row(str(metric), str(mode))
        console.print(metrics_table)

        defaults_table = self._add_columns(
            rich.table.Table(title="Default search space"),
            ["Hyperparameter", "Sample function", "Arguments", "Keyword arguments"],
        )
        for param, metadata in self._defaults.items():
            defaults_table.add_row(
                str(param),
                str(metadata["fn"]),
                str(metadata.get("args", [])),
                str(metadata.get("kwargs", {})),
            )
        console.print(defaults_table)

        if show_resources:
            resources = self._get_resources(available=True)
            resources_table = self._add_columns(
                rich.table.Table(title="Available resources"),
                ["Resource", "Quantity"],
            )
            resources_table.add_row("CPU cores", str(resources.get("CPU", "N/A")))
            resources_table.add_row("GPUs", str(resources.get("GPU", "N/A")))
            resources_table.add_row("Memory", str(resources.get("memory", "N/A")))
            console.print(resources_table)
