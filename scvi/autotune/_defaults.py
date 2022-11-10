from pytorch_lightning import LightningDataModule, LightningModule, Trainer

from scvi import model
from scvi.model.base import BaseModelClass
from scvi.module.base import BaseLatentModeModuleClass, BaseModuleClass
from scvi.train import TrainRunner

COLORS = [
    "dodger_blue1",
    "dark_violet",
    "green",
    "dark_orange",
]

TUNABLE_TYPE_TO_CLS = {
    "model": [BaseModelClass, BaseModuleClass, BaseLatentModeModuleClass],
    "train": [
        LightningDataModule,
        Trainer,
        TrainRunner,
    ],
    "train_plan": [
        LightningModule,
    ],
}

SUPPORTED = [model.SCVI]

DEFAULTS = {
    model.SCVI: {
        "lr": {"fn": "loguniform", "args": [1e-4, 1e-1]},
    }
}
