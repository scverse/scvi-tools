from .cellassign import CellAssign
from .gimvi import GIMVI
from .scar import SCAR
from .scbasset import SCBASSET
from .solo import SOLO
from .stereoscope import RNAStereoscope, SpatialStereoscope

# Jax modules
try:
    from .tangram import Tangram
except ImportError:
    class Tangram:
        def __init__(self, *args, **kwargs):
            raise NotImplementedError("This feature requires the 'jax' optional dependency.")

__all__ = [
    "SCAR",
    "SOLO",
    "GIMVI",
    "RNAStereoscope",
    "SpatialStereoscope",
    "CellAssign",
    "SCBASSET",
    "Tangram",
]
