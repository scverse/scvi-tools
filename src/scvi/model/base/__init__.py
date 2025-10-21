import warnings

from scvi import settings
from scvi.utils import error_on_missing_dependencies

from ._archesmixin import ArchesMixin
from ._base_model import (
    BaseMinifiedModeModelClass,
    BaseModelClass,
    BaseMudataMinifiedModeModelClass,
)
from ._differential import DifferentialComputation
from ._embedding_mixin import EmbeddingMixin
from ._pyromixin import (
    PyroJitGuideWarmup,
    PyroModelGuideWarmup,
    PyroSampleMixin,
    PyroSviTrainMixin,
)
from ._rnamixin import RNASeqMixin
from ._training_mixin import SemisupervisedTrainingMixin, UnsupervisedTrainingMixin
from ._vaemixin import VAEMixin

__all__ = [
    "ArchesMixin",
    "BaseModelClass",
    "RNASeqMixin",
    "VAEMixin",
    "UnsupervisedTrainingMixin",
    "SemisupervisedTrainingMixin",
    "PyroSviTrainMixin",
    "PyroSampleMixin",
    "PyroJitGuideWarmup",
    "PyroModelGuideWarmup",
    "DifferentialComputation",
    "BaseMinifiedModeModelClass",
    "BaseMudataMinifiedModeModelClass",
    "EmbeddingMixin",
]


def __getattr__(name: str):
    """Lazily provide object. If optional deps are missing, raise a helpful ImportError

    only when object is actually requested.
    """
    if name == "JaxTrainingMixin":
        warnings.warn(
            "In order to use the JaxTrainingMixin make sure to install scvi-tools[jax]",
            DeprecationWarning,
            stacklevel=settings.warnings_stacklevel,
        )

        error_on_missing_dependencies("flax", "jax", "jaxlib", "optax", "numpyro")
        from ._jaxmixin import JaxTrainingMixin as _JaxTrainingMixin

        return _JaxTrainingMixin
    raise AttributeError(f"module {__name__!r} has no attribute {name}")
