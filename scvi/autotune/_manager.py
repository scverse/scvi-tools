import inspect
import logging
import warnings
from collections import OrderedDict
from typing import Any, Callable, List, Literal, Optional, Tuple

import rich

from scvi._decorators import dependencies
from scvi._settings import settings
from scvi._types import AnnOrMuData
from scvi.autotune._defaults import (
    COLORS,
    DEFAULTS,
    SUPPORTED,
    TUNABLE_TYPE_TO_CLS,
    str_to_sample_fn,
)
from scvi.autotune._types import TunableMeta
from scvi.autotune._utils import in_notebook
from scvi.model.base import BaseModelClass

logger = logging.getLogger(__name__)


class TunerManager:
    """
    Provides an interface to validate and process a scvi-tools model class for autotuning.

    Validation of all inputs from :class:`~scvi.autotune.ModelTuner` is handled here.

    Parameters
    ----------
    model_cls
        :class:`~scvi.model.base.BaseModelClass` on which to tune hyperparameters. See
        :class:`~scvi.autotune.ModelTuner` for a list of supported models.
    """

    def __init__(self, model_cls: BaseModelClass):
        self._model_cls: BaseModelClass = self._validate_model(model_cls)
        self._registry: dict = self._validate_registry(model_cls)
        self._defaults: dict = self._validate_defaults(model_cls)

    @staticmethod
    def _validate_defaults(model_cls: BaseModelClass) -> dict:
        """Check whether the model class has a default search space defined."""
        if model_cls not in DEFAULTS:
            warnings.warn(
                f"No default search space available for {model_cls}.",
                UserWarning,
            )
        return DEFAULTS.get(model_cls, {})

    @staticmethod
    def _validate_model(model_cls: BaseModelClass) -> BaseModelClass:
        """Check whether the model class is supported."""
        if model_cls not in SUPPORTED:
            raise NotImplementedError(
                f"{model_cls} is currently unsupported, must be one of {SUPPORTED}."
            )
        return model_cls

    @staticmethod
    def _validate_registry(model_cls: BaseModelClass) -> dict:
        """
        Get tunable parameters and available metrics for a model class.

        For a given model class, checks whether a ``_tunables`` class property has been
        implemented. If so, iterates through the attribute and recursively discovers
        tunable parameters.

        Parameters
        ----------
        model_cls
            :class:`~scvi.model.base.BaseModelClass` on which to tune hyperparameters.

        Returns
        -------
        registry: dict
            Dictionary with the following keys:

            * ``'tunables'``: dictionary of tunable parameters
            * ``'metrics'``: dictionary of available metrics
        """

        def _cls_to_tunable_type(cls):
            for tunable_type, cls_list in TUNABLE_TYPE_TO_CLS.items():
                if any([issubclass(cls, c) for c in cls_list]):
                    return tunable_type

        def _get_tunables(
            attr: Any,
            parent: Any = None,
            tunable_type: Optional[str] = None,
        ) -> dict:
            tunables = {}
            if inspect.isfunction(attr):
                for k, v in inspect.signature(attr).parameters.items():
                    if isinstance(v.annotation, TunableMeta):
                        default = v.default
                        if default is inspect.Parameter.empty:
                            default = None
                        tunables[k] = {
                            "parent": parent,
                            "default": default,
                            "func": attr,
                            "tunable_type": tunable_type,
                        }
            elif inspect.isclass(attr) and hasattr(attr, "_tunables"):
                tunable_type = _cls_to_tunable_type(attr)
                for child in attr._tunables:
                    tunables.update(
                        _get_tunables(child, parent=attr, tunable_type=tunable_type)
                    )
            return tunables

        def _get_metrics(model_cls: BaseModelClass) -> dict:
            return {"validation_loss": "min"}

        registry = {
            "tunables": _get_tunables(model_cls),
            "metrics": _get_metrics(model_cls),
        }

        return registry

    def _parse_search_space(self, search_space: dict) -> Tuple[dict, dict]:
        """Parse a full search space configuration into separate dictionaries consumable by scvi-tools models."""
        model_kwargs = {}
        train_kwargs = {}
        plan_kwargs = {}
        _tunables = self._registry["tunables"]

        for k, v in search_space.items():
            _type = _tunables[k]["tunable_type"]
            if _type == "model":
                model_kwargs[k] = v
            elif _type == "train":
                train_kwargs[k] = v
            elif _type == "train_plan":
                plan_kwargs[k] = v
        train_kwargs["plan_kwargs"] = plan_kwargs

        return model_kwargs, train_kwargs

    def _validate_search_space(
        self,
        search_space: dict,
        use_defaults: bool,
        exclude: dict,
    ) -> dict:
        """Validate a user search space against the model class registry."""
        # validate user search space
        for key in search_space:
            if key not in self._registry["tunables"]:
                warnings.warn(
                    f"Provided parameter {key} is invalid for {self._model_cls}. Ignoring.",
                    UserWarning,
                )
                search_space.pop(key)

        # add defaults if requested
        if use_defaults:
            # e.g. convert "loguniform" to tune.loguniform
            defaults = {
                k: str_to_sample_fn(v["fn"])(*v["args"])
                for k, v in self._defaults.items()
            }
            logger.info(f"Initializing with default search space for {self._model_cls}")
            for param in exclude:
                if param not in defaults:
                    raise ValueError(
                        f"Excluded parameter {param} not in defaults. Ignoring."
                    )
                defaults.pop(param, None)
            search_space.update(defaults)

        return search_space

    def _validate_metrics(
        self,
        metric: str,
        additional_metrics: List[str],
    ) -> OrderedDict:
        """Validate a user metric(s) specification against the model class registry."""
        return OrderedDict({"validation_loss": "min"})

    @dependencies("ray.tune")
    def _validate_scheduler(
        self,
        scheduler: str,
        metrics: OrderedDict,
        scheduler_kwargs: dict,
    ) -> Any:
        """Validate a user scheduler specification and apply defaults."""
        from ray import tune

        metric = list(metrics.keys())[0]
        mode = metrics[metric]

        if scheduler == "asha":
            _default_kwargs = {
                "metric": metric,
                "mode": mode,
                "max_t": 100,
                "grace_period": 1,
                "reduction_factor": 2,
            }
            _scheduler = tune.schedulers.ASHAScheduler
        elif scheduler == "hyperband":
            _default_kwargs = {
                "metric": metric,
                "mode": mode,
                "max_t": 100,
                "reduction_factor": 2,
            }
            _scheduler = tune.schedulers.HyperBandScheduler
        elif scheduler == "median":
            _default_kwargs = {
                "metric": metric,
                "mode": mode,
                "grace_period": 1,
            }
            _scheduler = tune.schedulers.MedianStoppingRule
        elif scheduler == "pbt":
            _default_kwargs = {
                "metric": metric,
                "mode": mode,
            }
            _scheduler = tune.schedulers.PopulationBasedTraining
        elif scheduler == "fifo":
            _default_kwargs = {}
            _scheduler = tune.schedulers.FIFOScheduler

        _default_kwargs.update(scheduler_kwargs)
        return _scheduler(**_default_kwargs)

    @dependencies(["ray.tune", "hyperopt"])
    def _validate_searcher(
        self,
        searcher: str,
        metrics: OrderedDict,
        searcher_kwargs: dict,
    ) -> Any:
        """Validate a user searcher specification and apply defaults."""
        from ray import tune

        metric = list(metrics.keys())[0]
        mode = metrics[metric]

        if searcher in ["random", "grid"]:
            _default_kwargs = {
                "random_state": settings.seed,
            }
            _searcher = tune.search.basic_variant.BasicVariantGenerator
        elif searcher == "hyperopt":
            _default_kwargs = {
                "metric": metric,
                "mode": mode,
                "random_state_seed": settings.seed,
            }
            tune.search.SEARCH_ALG_IMPORT["hyperopt"]()  # tune not importing hyperopt
            _searcher = tune.search.hyperopt.HyperOptSearch

        _default_kwargs.update(searcher_kwargs)
        return _searcher(**_default_kwargs)

    def _validate_scheduler_and_searcher(
        self,
        scheduler: str,
        searcher: str,
        metrics: OrderedDict,
        scheduler_kwargs: dict,
        searcher_kwargs: dict,
    ) -> Tuple[Any, Any]:
        """Validate a user scheduler and searcher specifications for compatibility."""
        if scheduler not in ["asha", "hyperband", "median", "pbt", "fifo"]:
            raise ValueError(
                f"Provided scheduler {scheduler} is not supported. Must be "
                "one of ['asha', 'hyperband', 'median', 'pbt', 'fifo']."
            )
        if searcher not in ["random", "grid", "hyperopt"]:
            raise ValueError(
                f"Provided searcher {searcher} is not supported. Must be "
                "one of ['random', 'grid', 'hyperopt']."
            )
        if scheduler not in ["asha", "median", "hyperband"] and searcher not in [
            "random",
            "grid",
        ]:
            raise ValueError(
                f"Provided searcher {searcher} is incompatible with the "
                f"provided scheduler {scheduler}."
            )
        return (
            self._validate_scheduler(scheduler, metrics, scheduler_kwargs),
            self._validate_searcher(searcher, metrics, searcher_kwargs),
        )

    @dependencies("ray.tune")
    def _validate_reporter(
        self,
        reporter: bool,
        search_space: dict,
        metrics: OrderedDict,
    ) -> Any:
        """Validate a reporter depending on the user environment."""
        from ray import tune

        _metric_keys = list(metrics.keys())
        _search_keys = list(search_space.keys())
        _reporter_kwargs = {
            "metric_columns": _metric_keys,
            "parameter_columns": _search_keys,
            "metric": _metric_keys[0],
            "mode": metrics[_metric_keys[0]],
        }
        if not reporter:
            _reporter = None
        elif in_notebook():
            _reporter = tune.JupyterNotebookReporter(**_reporter_kwargs)
        else:
            _reporter = tune.CLIReporter(**_reporter_kwargs)
        return _reporter

    @dependencies("ray.tune")
    def _validate_trainable(
        self,
        adata: AnnOrMuData,
        metrics: OrderedDict,
        resources: dict,
        setup_kwargs: dict,
    ) -> Callable:
        """Create a trainable function consumable by :class:`~ray.tune.Tuner`."""
        from ray import tune
        from ray.tune.integration.pytorch_lightning import TuneReportCallback

        def _trainable(
            search_space: dict,
            model_cls: BaseModelClass = None,
            adata: AnnOrMuData = None,
            metric: str = None,
            setup_kwargs: dict = None,
        ) -> None:
            model_kwargs, train_kwargs = self._parse_search_space(search_space)
            model_cls.setup_anndata(adata, **setup_kwargs)
            model = model_cls(adata, **model_kwargs)
            monitor = TuneReportCallback(
                metric,
                on="validation_end",
            )
            model.train(
                check_val_every_n_epoch=1,
                callbacks=[monitor],
                enable_progress_bar=False,
                **train_kwargs,
            )

        return tune.with_resources(
            tune.with_parameters(
                _trainable,
                model_cls=self._model_cls,
                adata=adata,
                metric=list(metrics.keys())[0],
                setup_kwargs=setup_kwargs,
            ),
            resources=resources,
        )

    @dependencies(["ray.tune", "ray.air"])
    def get_tuner(
        self,
        adata: AnnOrMuData,
        *,
        metric: Optional[str] = None,
        additional_metrics: Optional[List[str]] = None,
        search_space: Optional[dict] = None,
        use_defaults: bool = True,
        exclude: Optional[dict] = None,
        scheduler: Literal["asha", "hyperband", "median", "pbt", "fifo"] = "asha",
        scheduler_kwargs: Optional[dict] = None,
        searcher: Literal["random", "grid", "hyperopt"] = "random",
        searcher_kwargs: Optional[dict] = None,
        reporter: bool = True,
        resources: Optional[dict] = None,
        setup_kwargs: Optional[dict] = None,
    ) -> Callable:
        """Configure a :class:`~ray.tune.Tuner` instance after validating user input."""
        from ray import air, tune

        additional_metrics = additional_metrics or []
        use_defaults = use_defaults if search_space is not None else True
        search_space = search_space or {}
        exclude = exclude or {}
        scheduler_kwargs = scheduler_kwargs or {}
        searcher_kwargs = searcher_kwargs or {}
        resources = resources or {}
        setup_kwargs = setup_kwargs or {}

        _search_space = self._validate_search_space(search_space, use_defaults, exclude)
        _metrics = self._validate_metrics(metric, additional_metrics)
        _scheduler, _searcher = self._validate_scheduler_and_searcher(
            scheduler,
            searcher,
            _metrics,
            scheduler_kwargs,
            searcher_kwargs,
        )
        _reporter = self._validate_reporter(reporter, _search_space, _metrics)
        _trainable = self._validate_trainable(
            adata,
            _metrics,
            resources,
            setup_kwargs,
        )
        tuner = tune.Tuner(
            trainable=_trainable,
            param_space=_search_space,
            tune_config=tune.TuneConfig(
                metric=list(_metrics.keys())[0],
                mode=list(_metrics.values())[0],
                scheduler=_scheduler,
                search_alg=_searcher,
                num_samples=5,
            ),
            run_config=air.RunConfig(
                name="scvi",
                progress_reporter=_reporter,
                failure_config=air.FailureConfig(
                    max_failures=2,
                ),
            ),
        )
        return tuner

    @staticmethod
    def _add_columns(table: rich.table.Table, columns: List[str], **kwargs):
        """Add columns to a :class:`~rich.table.Table` with default formatting."""
        _kwargs = {
            "justify": "center",
            "no_wrap": True,
            "overflow": "fold",
        }
        for i, column in enumerate(columns):
            table.add_column(column, style=COLORS[i], **_kwargs)
        return table

    def view_registry(self) -> None:
        """Prints a summary of the registry."""
        console = rich.console.Console(force_jupyter=in_notebook())

        # tunables
        tunables_table = self._add_columns(
            rich.table.Table(title="Tunable hyperparameters"),
            ["Hyperparameter", "Tunable type", "Default"],
        )
        for k, v in self._registry["tunables"].items():
            tunables_table.add_row(str(k), str(v["tunable_type"]), str(v["default"]))

        # metrics
        metrics_table = self._add_columns(
            rich.table.Table(title="Available metrics"), ["Metric", "Mode"]
        )
        for k, v in self._registry["metrics"].items():
            metrics_table.add_row(str(k), str(v))

        # defaults
        defaults_table = self._add_columns(
            rich.table.Table(title="Default search space"),
            ["Hyperparameter", "Sample type", "Arguments"],
        )
        for k, v in self._defaults.items():
            defaults_table.add_row(str(k), str(v["fn"]), str(v["args"]))

        # print
        console.print(f"Registry for {self._model_cls}")
        console.print(tunables_table)
        console.print(metrics_table)
        console.print(defaults_table)
