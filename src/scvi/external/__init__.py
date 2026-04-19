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
