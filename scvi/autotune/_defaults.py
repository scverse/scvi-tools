from pytorch_lightning import LightningDataModule, LightningModule, Trainer

from scvi import model
from scvi.model.base import BaseModelClass
from scvi.module.base import (
    BaseLatentModeModuleClass,
    BaseModuleClass,
    JaxBaseModuleClass,
)
from scvi.train import TrainRunner

COLORS = [
    "dodger_blue1",
    "dark_violet",
    "green",
    "dark_orange",
]

TUNABLE_TYPE_TO_CLS = {
    "model": [
        BaseModelClass,
        BaseModuleClass,
        BaseLatentModeModuleClass,
        JaxBaseModuleClass,
    ],
    "train": [
        LightningDataModule,
        Trainer,
        TrainRunner,
    ],
    "train_plan": [
        LightningModule,
    ],
}

SUPPORTED = [
    model.SCVI,
    model.TOTALVI,
    model.LinearSCVI,
    model.AUTOZI,
    model.SCANVI,
    model.PEAKVI,
    model.CondSCVI,
    model.DestVI,
    model.MULTIVI,
    model.AmortizedLDA,
    model.JaxSCVI,
]

DEFAULTS = {
    model.SCVI: {
        "lr": {"fn": "loguniform", "args": [1e-4, 1e-1]},
    },
    model.TOTALVI: {
        "lr": {"fn": "loguniform", "args": [1e-4, 1e-1]},
    },
    model.LinearSCVI: {
        "lr": {"fn": "loguniform", "args": [1e-4, 1e-1]},
    },
    model.AUTOZI: {
        "lr": {"fn": "loguniform", "args": [1e-4, 1e-1]},
    },
    model.SCANVI: {
        "lr": {"fn": "loguniform", "args": [1e-4, 1e-1]},
    },
    model.PEAKVI: {
        "lr": {"fn": "loguniform", "args": [1e-4, 1e-1]},
    },
    model.CondSCVI: {
        "lr": {"fn": "loguniform", "args": [1e-4, 1e-1]},
    },
    model.DestVI: {
        "lr": {"fn": "loguniform", "args": [1e-4, 1e-1]},
    },
    model.MULTIVI: {
        "lr": {"fn": "loguniform", "args": [1e-4, 1e-1]},
    },
    model.AmortizedLDA: {
        "lr": {"fn": "loguniform", "args": [1e-4, 1e-1]},
    },
    model.JaxSCVI: {
        "lr": {"fn": "loguniform", "args": [1e-4, 1e-1]},
    },
}
