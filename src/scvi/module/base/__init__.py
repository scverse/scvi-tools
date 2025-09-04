from scvi.utils import error_on_missing_dependencies

from ._base_module import (
    BaseMinifiedModeModuleClass,
    BaseModuleClass,
    LossOutput,
    PyroBaseModuleClass,
    SupervisedModuleClass,
)
from ._decorators import auto_move_data
from ._embedding_mixin import EmbeddingModuleMixin
from ._priors import GaussianPrior, MogPrior, VampPrior

__all__ = [
    "BaseModuleClass",
    "LossOutput",
    "PyroBaseModuleClass",
    "auto_move_data",
    "BaseMinifiedModeModuleClass",
    "SupervisedModuleClass",
    "EmbeddingModuleMixin",
    "GaussianPrior",
    "MogPrior",
    "VampPrior",
]


def __getattr__(name: str):
    """Lazily provide object. If optional deps are missing, raise a helpful ImportError

    only when object is actually requested.
    """
    if name == "flax_configure":
        error_on_missing_dependencies("flax", "jax", "jaxlib", "optax", "numpyro")
        from ._decorators import flax_configure as _flax_configure

        return _flax_configure
    if name == "JaxBaseModuleClass":
        error_on_missing_dependencies("flax", "jax", "jaxlib", "optax", "numpyro")
        from ._base_module import JaxBaseModuleClass as _JaxBaseModuleClass

        return _JaxBaseModuleClass
    if name == "TrainStateWithState":
        error_on_missing_dependencies("flax", "jax", "jaxlib", "optax", "numpyro")
        from ._base_module import TrainStateWithState as _TrainStateWithState

        return _TrainStateWithState
    raise AttributeError(f"module {__name__!r} has no attribute {name}")
