from collections import OrderedDict
from typing import Any, Callable, List, Optional, Tuple

import rich

from scvi._types import AnnOrMuData
from scvi.model.base import BaseModelClass


class TunerManager:
    """
    Internal manager for validation of inputs from :class:`~scvi.autotune.ModelTuner`.

    Parameters
    ----------
    model_cls
        :class:`~scvi.model.base.BaseModelClass` on which to tune
        hyperparameters. See :class:`~scvi.autotune.ModelTuner` for
        supported model classes.
    """

    def __init__(self, model_cls: BaseModelClass):
        pass

    def _validate_model(self, model_cls: BaseModelClass) -> None:
        """Checks if the model class is supported."""

    def _get_defaults(self, model_cls: BaseModelClass) -> dict:
        """Returns the model class's default search space if available."""

    def _get_registry(self, model_cls: BaseModelClass) -> dict:
        """
        Returns the model class's registry of tunable hyperparameters and metrics.

        For a given model class, checks whether a ``_tunables`` class property has been
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

    def _get_search_space(self, search_space: dict) -> Tuple[dict, dict]:
        """Parses a compact search space into separate kwargs dictionaries."""

    def _validate_search_space(
        self, search_space: dict, use_defaults: bool, exclude: List[str]
    ) -> dict:
        """Validates a search space against the hyperparameter registry."""

    def _validate_metrics(
        self, metric: str, additional_metrics: List[str]
    ) -> OrderedDict:
        """Validates a set of metrics against the metric registry."""

    def _validate_scheduler(
        self, scheduler: str, metrics: OrderedDict, scheduler_kwargs: dict
    ) -> Any:
        """Validates a trial scheduler."""

    def _validate_search_algorithm(
        self, searcher: str, metrics: OrderedDict, searcher_kwargs: dict
    ) -> Any:
        """Validates a hyperparameter search algorithm."""

    def _validate_scheduler_and_search_algorithm(
        self,
        scheduler: str,
        searcher: str,
        metrics: OrderedDict,
        scheduler_kwargs: dict,
        searcher_kwargs: dict,
    ) -> Tuple[Any, Any]:
        """Validates a scheduler and search algorithm pair for compatibility."""

    def _validate_reporter(
        self, reporter: bool, search_space: dict, metrics: OrderedDict
    ) -> Any:
        """Validates a reporter depending on the execution environment."""

    def _validate_resources(self, resources: dict) -> dict:
        """Validates a resources-use specification."""

    def _get_setup_kwargs(self, adata: AnnOrMuData) -> dict:
        """Retrieves the kwargs used for setting up ``adata`` with the model class."""

    def _get_trainable(
        self,
        adata: AnnOrMuData,
        metrics: OrderedDict,
        resources: dict,
        setup_kwargs: dict,
    ) -> Callable:
        """Returns a trainable function consumable by :class:`~ray.tune.Tuner`."""

    def _get_tuner(
        self,
        adata: AnnOrMuData,
        *,
        metric: Optional[str] = None,
        additional_metrics: Optional[List[str]] = None,
        search_space: Optional[dict] = None,
        num_samples: Optional[int] = None,
        use_defaults: bool = True,
        exclude: Optional[List[str]] = None,
        scheduler: Optional[str] = None,
        scheduler_kwargs: Optional[dict] = None,
        searcher: Optional[str] = None,
        searcher_kwargs: Optional[dict] = None,
        reporter: bool = True,
        resources: Optional[dict] = None,
    ) -> Any:
        """Configures a :class:`~ray.tune.Tuner` instance after validation."""

    def _add_columns(
        table: rich.table.Table, columns: List[str], **kwargs
    ) -> rich.table.Table:
        """Adds columns to a :class:`~rich.table.Table` with default formatting."""

    def _view_registry(self, shouw_resources: bool) -> None:
        """Displays a summary of the model class's registry and available resources."""
