from pytorch_lightning import LightningDataModule, LightningModule, Trainer

from scvi import model
from scvi.module.base import BaseModuleClass, JaxBaseModuleClass, PyroBaseModuleClass
from scvi.train import TrainRunner

# colors for rich table columns
COLORS = [
    "dodger_blue1",
    "dark_violet",
    "green",
    "dark_orange",
    "red",
]

# default rich table column kwargs
COLUMN_KWARGS = {
    "justify": "center",
    "no_wrap": True,
    "min_width": 10,
    "max_width": 50,
}

# maps classes to the type of hyperparameters they use
TUNABLE_TYPES = {
    "model": [
        BaseModuleClass,
        JaxBaseModuleClass,
        PyroBaseModuleClass,
    ],
    "train": [
        LightningDataModule,
        Trainer,
        TrainRunner,
    ],
    "training plan": [
        LightningModule,
    ],
}

# supported model classes
SUPPORTED = [
    model.SCVI,
    model.TOTALVI,
    model.SCANVI,
    model.MULTIVI,
    model.LinearSCVI,
    model.AUTOZI,
    model.DestVI,
    model.AmortizedLDA,
    model.PEAKVI,
    model.JaxSCVI,
]

# default hyperparameter search spaces for each model class
DEFAULTS = {
    model.SCVI: {
        "n_hidden": {"fn": "choice", "args": [[64, 128]]},
    }
}
