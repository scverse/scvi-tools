from scvi._decorators import dependencies

from ._types import Tunable, TunableMixin

__all__ = ["Tunable", "TunableMixin"]

try:
    from ._manager import TuneAnalysis, TunerManager
    from ._tuner import ModelTuner

    __all__ += ["ModelTuner", "TunerManager", "TuneAnalysis"]
except (ImportError, ModuleNotFoundError):

    @dependencies("ray.tune")
    def __getattr__(attr: str):
        pass
