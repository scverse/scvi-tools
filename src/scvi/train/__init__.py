from scvi.utils import error_on_missing_dependencies

from ._callbacks import (
    LoudEarlyStopping,
    SaveCheckpoint,
    ScibCallback,
)
from ._config import (
    AdversarialTrainingPlanConfig,
    ClassifierTrainingPlanConfig,
    KwargsConfig,
    LowLevelPyroTrainingPlanConfig,
    PyroTrainingPlanConfig,
    SemiSupervisedAdversarialTrainingPlanConfig,
    SemiSupervisedTrainingPlanConfig,
    TrainerConfig,
    TrainingPlanConfig,
    merge_kwargs,
)
from ._constants import METRIC_KEYS
from ._trainer import Trainer
from ._trainingplans import (
    AdversarialTrainingPlan,
    ClassifierTrainingPlan,
    LowLevelPyroTrainingPlan,
    PyroTrainingPlan,
    SemiSupervisedAdversarialTrainingPlan,
    SemiSupervisedTrainingPlan,
    TrainingPlan,
)
from ._trainrunner import TrainRunner

__all__ = [
    "merge_kwargs",
    "TrainingPlan",
    "TrainingPlanConfig",
    "Trainer",
    "TrainerConfig",
    "PyroTrainingPlan",
    "PyroTrainingPlanConfig",
    "LowLevelPyroTrainingPlan",
    "LowLevelPyroTrainingPlanConfig",
    "SemiSupervisedTrainingPlan",
    "SemiSupervisedTrainingPlanConfig",
    "SemiSupervisedAdversarialTrainingPlan",
    "SemiSupervisedAdversarialTrainingPlanConfig",
    "AdversarialTrainingPlan",
    "AdversarialTrainingPlanConfig",
    "ClassifierTrainingPlan",
    "ClassifierTrainingPlanConfig",
    "TrainRunner",
    "LoudEarlyStopping",
    "SaveCheckpoint",
    "ScibCallback",
    "METRIC_KEYS",
    "KwargsConfig",
]


def __getattr__(name: str):
    """Lazily provide object. If optional deps are missing, raise a helpful ImportError

    only when object is actually requested.
    """
    if name == "MlxTrainingPlan":
        error_on_missing_dependencies("mlx")
        from ._trainingplans import MlxTrainingPlan as _MlxTrainingPlan

        return _MlxTrainingPlan
    raise AttributeError(f"module {__name__!r} has no attribute {name}")
