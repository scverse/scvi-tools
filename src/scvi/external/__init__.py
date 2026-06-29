import warnings

from scvi import settings
from scvi.utils import error_on_missing_dependencies

from . import harreman
from .cellassign import CellAssign
from .contrastivevi import ContrastiveVI
from .cytovi import CYTOVI
from .decipher import Decipher
from .diagvi import DIAGVI
from .gimvi import GIMVI
from .harreman import HarremanAnalysis
from .joint_embedding_scvi import JointEmbeddingSCVI, JointEmbeddingVAE
from .methylvi import METHYLANVI, METHYLVI
from .mrvi import MRVI
from .poissonvi import POISSONVI
from .resolvi import RESOLVI
from .scar import SCAR
from .scbasset import SCBASSET
from .scviva import SCVIVA
from .solo import SOLO
from .stereoscope import RNAStereoscope, SpatialStereoscope
from .sysvi import SysVI
from .tangram import Tangram
from .totalanvi import TOTALANVI
from .velovi import VELOVI

__all__ = [
    "SCAR",
    "SOLO",
    "GIMVI",
    "JointEmbeddingSCVI",
    "JointEmbeddingVAE",
    "Decipher",
    "RNAStereoscope",
    "SpatialStereoscope",
    "CellAssign",
    "Tangram",
    "TOTALANVI",
    "SCBASSET",
    "POISSONVI",
    "ContrastiveVI",
    "SysVI",
    "VELOVI",
    "MRVI",
    "METHYLVI",
    "METHYLANVI",
    "RESOLVI",
    "SCVIVA",
    "CYTOVI",
    "DIAGVI",
    "Harreman",
    "HarremanAnalysis",
]
