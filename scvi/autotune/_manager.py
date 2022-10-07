import logging
import warnings
from inspect import signature

from scvi._types import AnnOrMuData
from scvi.model.base import BaseModelClass
from scvi.utils import attrdict

from ._defaults import SUPPORTED
from ._types import TunableMeta

logger = logging.getLogger(__name__)


class TunerManager:
    def __init__(self, model_cls: BaseModelClass, adata: AnnOrMuData):
        try:
            from ray import tune

            tune.choice([1])
        except ImportError:
            raise ImportError("Please install ray tune via `pip install 'ray[tune]'`.")

        self._model_cls: BaseModelClass = self._validate_model(model_cls)
        self._registry: attrdict = self._get_registry(model_cls, adata)

    def validate_search_space(self, search_config: dict, use_defaults) -> dict:
        """
        Given a user-provided search space, validate it against the registry.
        """
        logger.info(f"Using default parameters for {self._model_cls}")
        for key in search_config:
            if key not in self._registry:
                warnings.warn(
                    f"Search parameter {key} is invalid for the current model.",
                    UserWarning,
                )

    @staticmethod
    def _validate_model(model_cls: BaseModelClass) -> BaseModelClass:
        if model_cls not in SUPPORTED:
            raise NotImplementedError(
                f"ModelTuner currently unsupported for {model_cls}, must be one of "
                f"[{SUPPORTED}]."
            )
        return model_cls

    @staticmethod
    def _get_registry(model_cls: BaseModelClass, adata: AnnOrMuData) -> attrdict:
        """Get tunable parameters for a given model class."""
        dummy_model = model_cls(
            adata
        )  # TODO(martinkim0): If possible use @classproperty
        tunables = dict()
        for tunable_attr in dummy_model._tunables:
            for k, v in signature(tunable_attr.__init__).parameters.items():
                if isinstance(v.annotation, TunableMeta):
                    tunables[k] = tunable_attr

    def view_registry(self) -> None:
        self._view_registry(self._registry)

    @staticmethod
    def _view_registry(registry: dict) -> None:
        pass
