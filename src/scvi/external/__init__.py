from .cellassign import CellAssign
from .contrastivevi import ContrastiveVI
from .gimvi import GIMVI
<<<<<<< HEAD
from .methylvi import METHYLVI
=======
from .methylvi import MethylVI
>>>>>>> 28b58066 (Add documentation)
from .mrvi import MRVI
from .poissonvi import POISSONVI
from .scar import SCAR
from .scbasset import SCBASSET
from .solo import SOLO
from .stereoscope import RNAStereoscope, SpatialStereoscope
from .tangram import Tangram
from .velovi import VELOVI

__all__ = [
    "SCAR",
    "SOLO",
    "GIMVI",
    "RNAStereoscope",
    "SpatialStereoscope",
    "CellAssign",
    "Tangram",
    "SCBASSET",
    "POISSONVI",
    "ContrastiveVI",
    "VELOVI",
    "MRVI",
    "METHYLVI",
]
