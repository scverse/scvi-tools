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


def _invert_dict_with_lists(d: dict):
    d_inv = dict()
    for k, v in d.items():
        for v_ in v:
            d_inv[v_] = k


CLS_TO_TUNABLE_TYPE = _invert_dict_with_lists(TUNABLE_TYPE_TO_CLS)

SUPPORTED = [
    model.SCVI,
]

DEFAULTS = []
