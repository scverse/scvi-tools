import inspect
import logging
import sys
import warnings
from typing import Any, Optional

import rich

from scvi.autotune._defaults import DEFAULTS, SUPPORTED, TUNABLE_TYPE_TO_CLS
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
        See :class:`~scvi.autotune.ModelTuner` for supported models.
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

        def _cls_to_tunable_type(cls):
            for tunable_type, cls_list in TUNABLE_TYPE_TO_CLS.items():
                if any([issubclass(cls, c) for c in cls_list]):
                    return tunable_type
            return None

        def _get_tunables(
            attr: Any, parent: Optional[Any] = None, type: Optional[str] = None
        ) -> dict:
            tunables = dict()
            if inspect.isfunction(attr):
                for k, v in inspect.signature(attr).parameters.items():
                    if isinstance(v.annotation, TunableMeta):
                        default = v.default
                        if default is inspect.Parameter.empty:
                            default = None
                        tunables[k] = dict(
                            source=parent,
                            default=default,
                            func=attr,
                            type=type,
                        )
            elif inspect.isclass(attr) and hasattr(attr, "_tunables"):
                tunable_type = _cls_to_tunable_type(attr)
                for child in attr._tunables:
                    tunables.update(
                        _get_tunables(child, parent=attr, type=tunable_type)
                    )
            return tunables

        return attrdict(_get_tunables(model_cls), recursive=True)

    def validate_search_space(
        self, search_config: dict, use_defaults: bool, exclude: dict
    ) -> attrdict:
        """
        Given a user-provided search space, validate it against the registry.

        Parameters
        ----------
        search_config
            User-provided search space.
        use_defaults
            Whether to use default values for parameters not specified in the search.
        exclude
            If using defaults, whether to exclude certain parameters from the search.
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
        # Initialize search space with defaults if requested
        search_config_ = dict()
        if use_defaults:
            if model_cls in DEFAULTS:
                # TODO(martinkim0): Validate default params against registry
                logger.info(f"Initializing with default search space for {model_cls}")
                search_config_ = registry.copy()
            else:
                logger.info(f"No default search space available for {model_cls}.")

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

        # Separate model and train parameters
        model_config = {k: v for k, v in search_config_.items() if v.type == "model"}
        train_config = {k: v for k, v in search_config_.items() if v.type == "train"}

        return attrdict(dict(model=model_config, train=train_config), recursive=True)

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
        table.add_column(
            "Type",
            justify="center",
            style="yellow",
            no_wrap=True,
            overflow="fold",
        )
        for k, v in registry.items():
            table.add_row(str(k), str(v.source), str(v.default), str(v.type))

        console.print(table)
