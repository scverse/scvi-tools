from scvi.utils import error_on_missing_dependencies

from ._callbacks import (
    JaxModuleInit,
    LoudEarlyStopping,
    SaveBestState,
    SaveCheckpoint,
    SubSampleLabels,
)

__all__ = [
    "LoudEarlyStopping",
    "SaveBestState",
    "SaveCheckpoint",
    "JaxModuleInit",
    "SubSampleLabels",
]

try:
    error_on_missing_dependencies("scib_metrics")

    from ._scib_metrics import ScibCallback

    __all__.append("ScibCallback")
except ModuleNotFoundError:
    pass
