import inspect
import logging
import warnings
from collections import OrderedDict
from functools import partial
from typing import Any, Callable, List, Literal, Optional, Tuple

import rich

from scvi._decorators import dependencies
from scvi._settings import settings
from scvi._types import AnnOrMuData
from scvi.autotune._defaults import DEFAULTS, SUPPORTED, TUNABLE_TYPE_TO_CLS
from scvi.autotune._types import TunableMeta
from scvi.autotune._utils import in_notebook
from scvi.model.base import BaseModelClass

logger = logging.getLogger(__name__)


class TunerManager:
    """
    Provides an interface to validate and process a scvi-tools model class for use with
    :class:`~scvi.autotune.ModelTuner`.

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
        if model_cls not in DEFAULTS:
            warnings.warn(
                f"No default search space available for {model_cls}.",
                UserWarning,
            )
        return DEFAULTS.get(model_cls, {})

    @staticmethod
    def _validate_model(model_cls: BaseModelClass) -> BaseModelClass:
        """Validate input model class."""
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
            return None

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
            return {}

        registry = {
            "tunables": _get_tunables(model_cls),
            "metrics": _get_metrics(model_cls),
        }

        return registry

    def _parse_search_space(self, search_space: dict) -> Tuple[dict, dict]:
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

    @dependencies("ray.tune")
    def _validate_search_space(
        self, search_space: Optional[dict], use_defaults: bool, exclude: Optional[dict]
    ) -> dict:
        """Validate and process a user-provided search space."""
        # validate user search space
        search_space = search_space or {}
        exclude = exclude or {}
        for key in search_space:
            if key not in self._registry["tunables"]:
                warnings.warn(
                    f"Provided parameter {key} is invalid for {self._model_cls}. Ignoring.",
                    UserWarning,
                )
                search_space.pop(key)

        # add defaults if requested
        if use_defaults:
            defaults = self._defaults.copy()
            logger.info(f"Initializing with default search space for {self._model_cls}")
            for param in exclude:
                defaults.pop(param, None)
            search_space.update(defaults)

        return search_space

    def _validate_metrics(
        self,
        metric: Optional[str],
        additional_metrics: Optional[List[str]],
    ) -> OrderedDict:
        # metrics = additional_metrics or []
        # primary_metric = metric or None
        return OrderedDict({"validation_loss": "min"})

    def _validate_scheduler_and_searcher(
        self,
        scheduler: str,
        searcher: str,
        metrics: OrderedDict,
        scheduler_kwargs: Optional[dict],
        searcher_kwargs: Optional[dict],
    ) -> Tuple[Any, Any]:
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
            self._validate_scheduler(scheduler, metrics, scheduler_kwargs or {}),
            self._validate_searcher(searcher, metrics, searcher_kwargs or {}),
        )

    @dependencies("ray.tune")
    def _validate_scheduler(
        self, scheduler: str, metrics: OrderedDict, scheduler_kwargs: dict
    ) -> Any:
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
        self, searcher: str, metrics: OrderedDict, searcher_kwargs: dict
    ) -> None:
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

    @dependencies("ray.tune")
    def _validate_reporter(
        self, reporter: bool, search_space: dict, metrics: OrderedDict
    ) -> Any:
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
        from ray import tune
        from ray.tune.integration.pytorch_lightning import TuneReportCallback

        def _trainable(
            search_space: dict,
            model_cls: BaseModelClass = None,
            adata: AnnOrMuData = None,
            metric: str = None,
        ) -> None:
            model_cls.setup_anndata(adata, **setup_kwargs)
            model_kwargs, train_kwargs = self._parse_search_space(search_space)
            model = self._model_cls(adata, **model_kwargs)
            monitor = TuneReportCallback(
                metric,
                on="validation_end",
            )
            model.train(check_val_every_n_epoch=1, callbacks=[monitor] ** train_kwargs)

        return tune.with_resources(
            partial(
                _trainable,
                model_cls=self._model_cls,
                adata=adata,
                metric=list(metrics.keys())[0],
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
        """Configures a Ray Tuner instance."""
        from ray import air, tune

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
            resources or {},
            setup_kwargs or {},
        )
        print(_reporter)
        print(_searcher)
        print(_scheduler)
        print(_metrics)
        print(_search_space)
        tuner = tune.Tuner(
            trainable=_trainable,
            param_space=_search_space,
            tune_config=tune.TuneConfig(
                metric=list(_metrics.keys())[0],
                mode=list(_metrics.values())[0],
                scheduler=_scheduler,
                search_alg=_searcher,
                num_samples=None,
            ),
            run_config=air.RunConfig(
                name=None,
                progress_reporter=_reporter,
            ),
        )
        return tuner

    def view_registry(self) -> None:
        """Prints a summary of the registry."""
        console = rich.console.Console(force_jupyter=in_notebook())
        column_kwargs = {
            "justify": "center",
            "no_wrap": True,
            "overflow": "fold",
        }

        # tunables
        tunables_table = rich.table.Table(title="Tunable hyperparameters")
        columns = [
            ("Hyperparameter", "dodger_blue1"),
            ("Tunable type", "dark_violet"),
            ("Default", "green"),
        ]
        for column, color in columns:
            tunables_table.add_column(column, style=color, **column_kwargs)
        for k, v in self._registry["tunables"].items():
            tunables_table.add_row(str(k), str(v["tunable_type"]), str(v["default"]))

        # metrics
        metrics_table = rich.table.Table(title="Available metrics")
        columns = [
            ("Metric", "dodger_blue1"),
            ("Mode", "dark_violet"),
        ]
        for column, color in columns:
            metrics_table.add_column(column, style=color, **column_kwargs)
        metrics = self._registry["metrics"]
        for k, v in metrics.items():
            metrics_table.add_row(str(k), str(v))

        # defaults
        defaults_table = rich.table.Table(title="Default search space")
        columns = [
            ("Hyperparameter", "dodger_blue1"),
            ("Sample type", "dark_violet"),
            ("Arguments", "green"),
        ]
        for column, color in columns:
            defaults_table.add_column(column, style=color, **column_kwargs)
        for k, v in self._defaults.items():
            defaults_table.add_row(str(k), str(v["fn"]), str(v["args"]))

        # print
        console.print(f"Registry for {self._model_cls}")
        console.print(tunables_table)
        console.print(metrics_table)
        console.print(defaults_table)
