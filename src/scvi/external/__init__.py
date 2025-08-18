from scvi.utils import error_on_missing_dependencies

from .cellassign import CellAssign
from .contrastivevi import ContrastiveVI
from .decipher import Decipher
from .gimvi import GIMVI
from .methylvi import METHYLANVI, METHYLVI
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
    "METHYLVI",
    "METHYLANVI",
    "RESOLVI",
    "SCVIVA",
]


def __getattr__(name: str):
    """
    Lazily provide object. If optional deps are missing, raise a helpful ImportError

    only when object is actually requested.
    """
    if name == "MRVI":
        error_on_missing_dependencies("flax", "jax", "jaxlib", "optax", "numpyro", "xarray")
        from .mrvi_torch import MRVI as _MRVI

        return _MRVI
    if name == "TorchMRVI":
        # TODO: REMOVE THIS ONCE TORCHMRVI IS READY AND MOVE TO THE __all__ above
        error_on_missing_dependencies("flax", "jax", "jaxlib", "optax", "numpyro", "xarray")
        from .mrvi_torch import TorchMRVI as _TorchMRVI

        return _TorchMRVI
    if name == "Tangram":
        error_on_missing_dependencies("flax", "jax", "jaxlib", "optax", "numpyro")
        from .tangram import Tangram as _Tangram

        return _Tangram
    raise AttributeError(f"module {__name__!r} has no attribute {name}")
