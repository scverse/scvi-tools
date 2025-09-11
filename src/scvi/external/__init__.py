import warnings

from scvi import settings
from scvi.utils import error_on_missing_dependencies

from .cellassign import CellAssign
from .contrastivevi import ContrastiveVI
from .cytovi import CYTOVI
from .decipher import Decipher
from .gimvi import GIMVI
from .methylvi import METHYLANVI, METHYLVI
from .mrvi import MRVI
from .mrvi_torch import TorchMRVI
from .poissonvi import POISSONVI
from .resolvi import RESOLVI
from .scar import SCAR
from .scbasset import SCBASSET
from .scviva import SCVIVA
from .solo import SOLO
from .stereoscope import RNAStereoscope, SpatialStereoscope
from .sysvi import SysVI
from .totalanvi import TOTALANVI
from .velovi import VELOVI

__all__ = [
    "SCAR",
    "SOLO",
    "GIMVI",
    "Decipher",
    "RNAStereoscope",
    "SpatialStereoscope",
    "CellAssign",
    "TOTALANVI",
    "SCBASSET",
    "POISSONVI",
    "ContrastiveVI",
    "SysVI",
    "VELOVI",
    "MRVI",
    "TorchMRVI",
    "METHYLVI",
    "METHYLANVI",
    "RESOLVI",
    "SCVIVA",
    "CYTOVI",
]


def __getattr__(name: str):
    """
    Lazily provide object. If optional deps are missing, raise a helpful ImportError

    only when object is actually requested.
    """
    if name == "JaxMRVI":
        warnings.warn(
            "In order to use the Jax version of MRVI make sure to install scvi-tools[jax]",
            DeprecationWarning,
            stacklevel=settings.warnings_stacklevel,
        )

        error_on_missing_dependencies("flax", "jax", "jaxlib", "optax", "numpyro")
        from .mrvi_jax import JaxMRVI as _JaxMRVI

        return _JaxMRVI

    if name == "Tangram":
        warnings.warn(
            "In order to use the TANGRAM make sure to install scvi-tools[jax]",
            DeprecationWarning,
            stacklevel=settings.warnings_stacklevel,
        )

        error_on_missing_dependencies("flax", "jax", "jaxlib", "optax", "numpyro")
        from .tangram import Tangram as _Tangram

        return _Tangram
    raise AttributeError(f"module {__name__!r} has no attribute {name}")
