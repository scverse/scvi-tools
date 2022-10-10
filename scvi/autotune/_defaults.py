from pytorch_lightning import LightningDataModule, LightningModule, Trainer

from scvi import model
from scvi.model.base import BaseModelClass
from scvi.module.base import BaseLatentModeModuleClass, BaseModuleClass
from scvi.train import TrainRunner

TUNABLE_TYPE_TO_CLS = dict(
    model=[BaseModelClass, BaseModuleClass, BaseLatentModeModuleClass],
    train=[
        LightningDataModule,
        Trainer,
        TrainRunner,
    ],
    train_plan=[
        LightningModule,
    ],
)

SUPPORTED = [
    model.SCVI,
]

DEFAULTS = {model.SCVI: dict()}
