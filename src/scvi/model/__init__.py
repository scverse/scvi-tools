from scvi.utils import error_on_missing_dependencies

from . import utils
from ._amortizedlda import AmortizedLDA
from ._autozi import AUTOZI
from ._condscvi import CondSCVI
from ._destvi import DestVI
from ._linear_scvi import LinearSCVI
from ._multivi import MULTIVI
from ._peakvi import PEAKVI
from ._scanvi import SCANVI
from ._scvi import SCVI
from ._totalvi import TOTALVI

__all__ = [
    "SCVI",
    "TOTALVI",
    "LinearSCVI",
    "AUTOZI",
    "SCANVI",
    "PEAKVI",
    "CondSCVI",
    "DestVI",
    "MULTIVI",
    "AmortizedLDA",
    "utils",
]


def __getattr__(name: str):
    """Lazily provide object. If optional deps are missing, raise a helpful ImportError

    only when object is actually requested.
    """
    if name == "JaxSCVI":
        error_on_missing_dependencies("flax", "jax", "jaxlib", "optax", "numpyro", "xarray")
        from ._jaxscvi import JaxSCVI as _JaxSCVI

        return _JaxSCVI
    raise AttributeError(f"module {__name__!r} has no attribute {name}")
