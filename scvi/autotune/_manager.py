import inspect
import logging
import warnings
from collections import OrderedDict
from typing import Any, Callable, List, Optional, Tuple

import rich

try:
    from ray import air, tune
    from ray.tune.integration.pytorch_lightning import TuneReportCallback
except ImportError:
    pass

from scvi._decorators import dependencies
from scvi._types import AnnOrMuData
from scvi.data._constants import _SETUP_ARGS_KEY, _SETUP_METHOD_NAME
from scvi.model.base import BaseModelClass

from ._defaults import COLORS, COLUMN_KWARGS, DEFAULTS, SUPPORTED, TUNABLE_TYPES
from ._types import TunableMeta
from ._utils import in_notebook

logger = logging.getLogger(__name__)


class TunerManager:
    """
    Internal manager for validation of inputs from :class:`~scvi.autotune.ModelTuner`.

    Parameters
    ----------
    model_cls
        :class:`~scvi.model.base.BaseModelClass` on which to tune hyperparameters. See
        :class:`~scvi.autotune.ModelTuner` forsupported model classes.
    """

    def __init__(self, model_cls: BaseModelClass):
        self._model_cls: BaseModelClass = self._validate_model_cls(model_cls)
        self._defaults: dict = self._get_defaults(self._model_cls)
        self._registry: dict = self._get_registry(self._model_cls)

    def _validate_model_cls(self, model_cls: BaseModelClass) -> BaseModelClass:
        """Checks if the model class is suppo rted."""
        if model_cls not in SUPPORTED:
            raise NotImplementedError(
                f"{model_cls} is currently unsupported. Please see ModelTuner for a "
                "list of supported model classes."
            )
        return model_cls

    def _get_defaults(self, model_cls: BaseModelClass) -> dict:
        """Returns the model class's default search space if available."""
        if model_cls not in DEFAULTS:
            warnings.warn(
                f"No default search space available for {model_cls}.",
                UserWarning,
            )
        return DEFAULTS.get(model_cls, {})

    def _get_registry(self, model_cls: BaseModelClass) -> dict:
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
            return "unknown"

        def _get_tunables(
            attr: Any, parent: Any = None, tunable_type: Optional[str] = None
        ) -> dict:
            tunables = {}
            if inspect.isfunction(attr):
                # check if function kwargs are tunable
                for kwarg, metadata in inspect.signature(attr).parameters.items():
                    if not isinstance(metadata.annotation, TunableMeta):
                        continue
                    default_val = metadata.default
                    if default_val is inspect.Parameter.empty:
                        default_val = None
                    tunables[kwarg] = {
                        "parent_class": parent,
                        "default_value": default_val,
                        "function": attr,
                        "tunable_type": tunable_type,
                    }
            elif inspect.isclass(attr) and hasattr(attr, "_tunables"):
                # recursively check if `_tunables` is implemented
                tunable_type = _cls_to_tunable_type(attr)
                for child in attr._tunables:
                    tunables.update(
                        _get_tunables(child, parent=attr, tunable_type=tunable_type)
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
    def _validate_search_space(
        self, search_space: dict, use_defaults: bool, exclude: List[str]
    ) -> dict:
        """Validates a search space against the hyperparameter registry."""
        # validate user-provided search space
        for param in search_space:
            if param in self._registry["tunables"]:
                continue
            warnings.warn(
                f"Provided parameter {param} is invalid for {self._model_cls.__name__}."
                " Please see available parameters with `ModelTuner.info()`. "
                "Ignoring parameter.",
                UserWarning,
            )
            search_space.pop(param)

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
            for param in exclude:
                if param not in _search_space:
                    warnings.warn(
                        f"Excluded parameter {param} not in defaults search space. "
                        "Ignoring parameter.",
                        UserWarning,
                    )
                _search_space.pop(param, None)

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
                f"Provided metric {metric} is invalid for {self._model_cls.__name__}. "
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

    @dependencies("ray.tune")
    def _validate_scheduler(
        self, scheduler: str, metrics: OrderedDict, scheduler_kwargs: dict
    ) -> Any:
        """Validates a trial scheduler."""
        metric = list(metrics.keys())[0]
        mode = metrics[metric]
        _kwargs = {"metric": metric, "mode": mode}

        if scheduler == "asha":
            _default_kwargs = {
                "max_t": 100,
                "grace_period": 1,
                "reduction_factor": 2,
            }
            _scheduler = tune.schedulers.ASHAScheduler
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
        metric = list(metrics.keys())[0]
        mode = metrics[metric]

        if searcher == "random":
            _default_kwargs = {}
            _searcher = tune.search.basic_variant.BasicVariantGenerator
        elif searcher == "grid":
            _default_kwargs = {}
            _searcher = tune.search.basic_variant.BasicVariantGenerator
        elif searcher == "hyperopt":
            _default_kwargs = {
                "metric": metric,
                "mode": mode,
            }
            tune.search.SEARCH_ALG_IMPORT["hyperopt"]()  # tune not importing hyperopt
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
        if scheduler not in ["asha", "hyperband", "median", "pbt", "fifo"]:
            raise ValueError(
                f"Provided scheduler {scheduler} is unsupported. Must be one of  "
                "['asha', 'hyperband', 'median', 'pbt', 'fifo']. ",
            )
        if searcher not in ["random", "grid", "hyperopt"]:
            raise ValueError(
                f"Provided searcher {searcher} is unsupported. Must be one of "
                "['random', 'grid', 'hyperopt']. ",
            )
        if scheduler not in ["asha", "median", "hyperband"] and searcher not in [
            "random",
            "grid",
        ]:
            raise ValueError(
                f"Provided scheduler {scheduler} is incompatible with the provided "
                f"searcher {searcher}. Please see "
                "https://docs.ray.io/en/latest/tune/key-concepts.html for more info."
            )

        _scheduler = self._validate_scheduler(scheduler, metrics, scheduler_kwargs)
        _searcher = self._validate_search_algorithm(searcher, metrics, searcher_kwargs)
        return _scheduler, _searcher

    @dependencies("ray.tune")
    def _validate_reporter(
        self, reporter: bool, search_space: dict, metrics: OrderedDict
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
        ) -> None:
            model_kwargs, train_kwargs = self._get_search_space(search_space)
            # TODO: generalize to models with mudata
            getattr(model_cls, setup_method_name)(adata, **setup_kwargs)
            model = model_cls(adata, **model_kwargs)
            monitor = TuneReportCallback(metric, on="validation_end")
            # TODO: adaptive max_epochs
            model.train(
                max_epochs=max_epochs,
                check_val_every_n_epoch=1,
                callbacks=[monitor],
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
        )
        return tune.with_resources(_wrap_params, resources=resources)

    @dependencies(["ray.tune", "ray.air"])
    def _get_tuner(
        self,
        adata: AnnOrMuData,
        *,
        metric: Optional[str] = None,
        additional_metrics: Optional[List[str]] = None,
        search_space: Optional[dict] = None,
        use_defaults: bool = True,
        exclude: Optional[List[str]] = None,
        num_samples: Optional[int] = None,
        max_epochs: Optional[int] = None,
        scheduler: Optional[str] = None,
        scheduler_kwargs: Optional[dict] = None,
        searcher: Optional[str] = None,
        searcher_kwargs: Optional[dict] = None,
        reporter: bool = True,
        resources: Optional[dict] = None,
    ) -> Any:
        """Configures a :class:`~ray.tune.Tuner` instance after validation."""
        metric = metric or list(self._registry["metrics"].keys())[0]
        additional_metrics = additional_metrics or []
        search_space = search_space or {}
        exclude = exclude or []
        num_samples = num_samples or 10
        max_epochs = max_epochs or 10
        scheduler = scheduler or "asha"
        scheduler_kwargs = scheduler_kwargs or {}
        searcher = searcher or "hyperopt"
        searcher_kwargs = searcher_kwargs or {}
        resources = resources or {}

        _ = self._model_cls(adata)
        _metrics = self._validate_metrics(metric, additional_metrics)
        _search_space = self._validate_search_space(search_space, use_defaults, exclude)
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

        tune_config = tune.tune_config.TuneConfig(
            scheduler=_scheduler,
            search_alg=_searcher,
            num_samples=num_samples,
        )
        # TODO: add kwarg for name or auto-generate name?
        run_config = air.config.RunConfig(
            name="scvi-tune",
            progress_reporter=_reporter,
        )
        tuner = tune.Tuner(
            trainable=_trainable,
            param_space=_search_space,
            tune_config=tune_config,
            run_config=run_config,
        )
        return tuner

    def _add_columns(
        self, table: rich.table.Table, columns: List[str]
    ) -> rich.table.Table:
        """Adds columns to a :class:`~rich.table.Table` with default formatting."""
        for i, column in enumerate(columns):
            table.add_column(column, style=COLORS[i], **COLUMN_KWARGS)
        return table

    def _view_registry(self, show_resources: bool) -> None:
        """Displays a summary of the model class's registry and available resources."""
        console = rich.console.Console(force_jupyter=in_notebook())

        tunables_table = self._add_columns(
            rich.table.Table(title="Tunable hyperparameters"),
            ["Hyperparameter", "Tunable type", "Default value", "Source"],
        )
        for param, metadata in self._registry["tunables"].items():
            tunables_table.add_row(
                str(param),
                str(metadata["tunable_type"]),
                str(metadata["default_value"]),
                str(metadata["parent_class"]),
            )

        metrics_table = self._add_columns(
            rich.table.Table(title="Available metrics"),
            ["Metric", "Mode"],
        )
        for metric, mode in self._registry["metrics"].items():
            metrics_table.add_row(str(metric), str(mode))

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

        console.print(f"Registry for {self._model_cls.__name__}")
        console.print(tunables_table)
        console.print(metrics_table)
        console.print(defaults_table)

        if show_resources:
            # TODO: retrieve available resources
            pass
