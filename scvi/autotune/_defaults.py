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
]

# default rich table column kwargs
COLUMN_KWARGS = {
    "justify": "center",
    "no_wrap": True,
    "overflow": "fold",
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
    "train_plan": [
        LightningModule,
    ],
}

# supported model classes
SUPPORTED = [model.SCVI]

# default hyperparameter search spaces for each model class
DEFAULTS = {
    model.SCVI: {
        "n_hidden": {"fn": "choice", "args": [[64, 128]]},
    }
}
