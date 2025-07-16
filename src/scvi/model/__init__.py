from scvi.utils import is_package_installed

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

if is_package_installed("numpyro") and is_package_installed("jax"):
    from ._jaxscvi import JaxSCVI

    __all__ += ["JaxSCVI"]
