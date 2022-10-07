import inspect
import logging
import sys
import warnings
from typing import Any

import rich

from scvi.autotune._defaults import SUPPORTED
from scvi.autotune._types import TunableMeta
from scvi.model.base import BaseModelClass
from scvi.utils import attrdict

logger = logging.getLogger(__name__)


class TunerManager:
    """
    Provides an interface to validate and process a scvi-tools model class for use with
    :class:`~scvi.autotune.ModelTuner`.

    Validation of all inputs from :class:`~scvi.autotune.ModelTuner` methods is handled
    in this class.

    Parameters
    ----------
    model_cls
        :class:`~scvi.model.base.BaseModelClass` on which to tune hyperparameters.
        Currently supports one of the following:

        * :class:`~scvi.model.SCVI`
    """

    def __init__(self, model_cls: BaseModelClass):
        try:
            from ray import tune

            tune.choice([None])  # for static code check
        except ImportError:
            raise ImportError("Please install ray tune via `pip install 'ray[tune]'`.")

        self._model_cls = self._validate_model(model_cls)
        self._registry = self._get_registry(model_cls)

    @staticmethod
    def _validate_model(model_cls: BaseModelClass) -> BaseModelClass:
        if model_cls not in SUPPORTED:
            raise NotImplementedError(
                f"ModelTuner currently unsupported for {model_cls}, must be one of "
                f"[{SUPPORTED}]."
            )
        return model_cls

    @staticmethod
    def _get_registry(model_cls: BaseModelClass) -> attrdict:
        """
        Get tunable parameters for a given model class.

        For a given model class, check whether a ``_tunables`` class property has been
        implemented. If so, iterate through the attribute and recursively find tunable
        parameters.

        Parameters
        ----------
        model_cls
            Model class to check for tunable parameters.

        Returns
        -------
        tunables: attrdict
            Dictionary of tunable parameters, where keys are parameter names and values
            are the sources of the parameters.
        """

        def _get_tunables(attr: Any) -> dict:
            tunables = dict()
            if hasattr(attr, "_tunables"):
                for child in attr._tunables:
                    if inspect.isfunction(child):
                        for k, v in inspect.signature(attr).parameters.items():
                            if isinstance(v.annotation, TunableMeta):
                                default = v.default
                                if default is inspect.Parameter.empty:
                                    default = None
                                tunables[k] = (attr, default, child)
                    else:
                        tunables.update(_get_tunables(child))
            return tunables

        return attrdict(_get_tunables(model_cls))

    def validate_search_space(
        self, search_config: dict, use_defaults: bool, exclude: dict
    ) -> attrdict:
        """
        Given a user-provided search space, validate it against the registry.
        """
        return self._validate_search_space(
            self._model_cls,
            self._registry,
            search_config,
            use_defaults,
            exclude,
        )

    @staticmethod
    def _validate_search_space(
        model_cls: BaseModelClass,
        registry: attrdict,
        search_config: dict,
        use_defaults: bool,
        exclude: dict,
    ) -> attrdict:
        if use_defaults:
            logger.info(f"Incorporating default search space for {model_cls}")
            search_config_ = registry.copy()
        else:
            search_config_ = dict()

        # Validate user search space
        for key in search_config:
            if key not in registry:
                warnings.warn(
                    f"Provided parameter {key} is invalid for {model_cls}. Ignoring.",
                    UserWarning,
                )
                _ = search_config.pop(key, None)
        search_config_.update(search_config)

        # Validate excluded parameters
        for key in exclude:
            if key not in search_config_:
                warnings.warn(
                    f"Excluded parameter {key} is invalid for {model_cls}. Ignoring.",
                )
            _ = search_config_.pop(key, None)

        return attrdict(search_config_)

    def view_registry(self) -> None:
        """
        Prints a summary of the registry.
        """
        self._view_registry(self._registry, self._model_cls)

    @staticmethod
    def _view_registry(registry: attrdict, model_cls: BaseModelClass) -> None:
        in_colab = "google.colab" in sys.modules
        force_jupyter = None if not in_colab else True
        console = rich.console.Console(force_jupyter=force_jupyter)

        table = rich.table.Table(title=f"ModelTuner registry for {model_cls}")
        table.add_column(
            "Hyperparameter name",
            justify="center",
            style="dodger_blue1",
            no_wrap=True,
            overflow="fold",
        )
        table.add_column(
            "Source",
            justify="center",
            style="dark_violet",
            no_wrap=True,
            overflow="fold",
        )
        table.add_column(
            "Default value",
            justify="center",
            style="green",
            no_wrap=True,
            overflow="fold",
        )
        for k, v in registry.items():
            table.add_row(str(k), str(v[0]), str(v[1]))

        console.print(table)
