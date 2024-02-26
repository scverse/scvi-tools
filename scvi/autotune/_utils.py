import inspect
from typing import Any, Literal

from lightning.pytorch import LightningDataModule, LightningModule, Trainer
from rich.table import Table

from scvi.module.base import BaseModuleClass, JaxBaseModuleClass, PyroBaseModuleClass
from scvi.train import TrainRunner

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
    "training_plan": [
        LightningModule,
    ],
}

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


def _cls_to_tunable_type(cls: type) -> Literal["model", "train", "training_plan"]:
    """Get the tunable type of a class."""
    for tunable_type, cls_list in TUNABLE_TYPES.items():
        if any(issubclass(cls, c) for c in cls_list):
            return tunable_type

    return None


def _get_tunable_function_params(fn: callable, parent: type, tunable_type: str) -> dict:
    from scvi._types import TunableMeta

    tunable_params = {}
    for param, metadata in inspect.signature(fn).parameters.items():
        cond1 = isinstance(metadata.annotation, TunableMeta)
        cond2 = "Tunable" in str(metadata.annotation)  # postponed annotation evaluation
        if not cond1 and not cond2:
            continue

        default_val = None
        if metadata.default != inspect.Parameter.empty:
            default_val = metadata.default

        if cond1:
            annotation = metadata.annotation.__args__[0]
            if hasattr(annotation, "__args__"):
                # e.g. if type is Literal, get its arguments
                annotation = annotation.__args__
            elif hasattr(annotation, "__name__"):
                annotation = annotation.__name__
            else:
                annotation = str(annotation)
        else:
            annotation = metadata.annotation
            annotation = annotation[annotation.find("[") + 1 : annotation.rfind("]")]

        tunable_params[param] = {
            "tunable_type": tunable_type,
            "default_value": default_val,
            "source": parent.__name__,
            "annotation": annotation,
        }

    return tunable_params


def get_tunables(attr: Any, parent: Any = None, tunable_type: str | None = None) -> dict:
    tunables = {}
    if inspect.isfunction(attr):
        return _get_tunable_function_params(attr, parent, tunable_type)
    for child in getattr(attr, "_tunables", {}):
        tunables.update(get_tunables(child, parent=attr, tunable_type=_cls_to_tunable_type(attr)))
    return tunables


def add_rich_table_columns(table: Table, columns: list[str]) -> Table:
    """Add columns to a :class:`~rich.table.Table using default styles."""
    for i, column in enumerate(columns):
        table.add_column(column, style=COLORS[i % len(COLORS)], **COLUMN_KWARGS)
    return table
