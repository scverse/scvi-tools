import warnings

from scvi import settings

try:
    from ray import tune
except ModuleNotFoundError as e:
    raise ModuleNotFoundError(
        "Please install ray[tune] to use scvi.autotune. Skipping import."
    ) from e


from ._manager import TuneAnalysis, TunerManager
from ._tuner import ModelTuner

__all__ = ["ModelTuner", "TunerManager", "TuneAnalysis"]
