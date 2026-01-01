from __future__ import annotations

from collections.abc import Callable, Iterable, Mapping
from dataclasses import dataclass, field
from typing import Any, Protocol, TYPE_CHECKING, runtime_checkable

import torch

from scvi._constants import REGISTRY_KEYS

if TYPE_CHECKING:
    from typing import Literal

    from lightning.pytorch.accelerators import Accelerator
    from lightning.pytorch.loggers import Logger

type TorchOptimizerCreator = Callable[[Iterable[torch.Tensor]], torch.optim.Optimizer]


@runtime_checkable
class KwargsConfig(Protocol):
    """Protocol for config objects that can be expanded into kwargs."""

    def to_kwargs(self) -> dict[str, Any]:
        """Return keyword arguments compatible with a downstream constructor."""


type KwargsLike = Mapping[str, Any] | KwargsConfig


def _coerce_kwargs(value: KwargsLike | None, *, name: str) -> dict[str, Any]:
    """Normalize a kwargs-like object into a dict."""
    if value is None:
        return {}
    if isinstance(value, dict):
        return value
    if isinstance(value, Mapping):
        return dict(value)
    to_kwargs = getattr(value, "to_kwargs", None)
    if callable(to_kwargs):
        out = to_kwargs()
        if not isinstance(out, dict):
            raise TypeError(f"{name}.to_kwargs() must return a dict, got {type(out)!r}.")
        return out
    raise TypeError(f"{name} must be a mapping or a config with to_kwargs().")


def merge_kwargs(
    config: KwargsLike | None,
    overrides: KwargsLike | None,
    *,
    name: str,
) -> dict[str, Any]:
    """Merge config kwargs with overrides, with overrides taking precedence."""
    merged = dict(_coerce_kwargs(config, name=f"{name}_config"))
    merged.update(_coerce_kwargs(overrides, name=name))
    return merged


@dataclass
class TrainingPlanConfig:
    """Config for :class:`~scvi.train.TrainingPlan`."""

    optimizer: Literal["Adam", "AdamW", "Custom"] = "Adam"
    optimizer_creator: TorchOptimizerCreator | None = None
    lr: float = 1e-3
    update_only_decoder: bool = False
    weight_decay: float = 1e-6
    eps: float = 0.01
    n_steps_kl_warmup: int | None = None
    n_epochs_kl_warmup: int | None = 400
    reduce_lr_on_plateau: bool = False
    lr_factor: float = 0.6
    lr_patience: int = 30
    lr_threshold: float = 0.0
    lr_scheduler_metric: Literal[
        "elbo_validation",
        "reconstruction_loss_validation",
        "kl_local_validation",
    ] = "elbo_validation"
    lr_min: float = 0.0
    max_kl_weight: float = 1.0
    min_kl_weight: float = 0.0
    compile: bool = False
    compile_kwargs: dict | None = None
    on_step: bool | None = False
    on_epoch: bool | None = True
    loss_kwargs: dict[str, Any] = field(default_factory=dict)

    def to_kwargs(self) -> dict[str, Any]:
        kwargs = {
            "optimizer": self.optimizer,
            "optimizer_creator": self.optimizer_creator,
            "lr": self.lr,
            "update_only_decoder": self.update_only_decoder,
            "weight_decay": self.weight_decay,
            "eps": self.eps,
            "n_steps_kl_warmup": self.n_steps_kl_warmup,
            "n_epochs_kl_warmup": self.n_epochs_kl_warmup,
            "reduce_lr_on_plateau": self.reduce_lr_on_plateau,
            "lr_factor": self.lr_factor,
            "lr_patience": self.lr_patience,
            "lr_threshold": self.lr_threshold,
            "lr_scheduler_metric": self.lr_scheduler_metric,
            "lr_min": self.lr_min,
            "max_kl_weight": self.max_kl_weight,
            "min_kl_weight": self.min_kl_weight,
            "compile": self.compile,
            "compile_kwargs": self.compile_kwargs,
            "on_step": self.on_step,
            "on_epoch": self.on_epoch,
        }
        kwargs.update(self.loss_kwargs)
        return kwargs


@dataclass
class AdversarialTrainingPlanConfig:
    """Config for :class:`~scvi.train.AdversarialTrainingPlan`."""

    optimizer: Literal["Adam", "AdamW", "Custom"] = "Adam"
    optimizer_creator: TorchOptimizerCreator | None = None
    lr: float = 1e-3
    weight_decay: float = 1e-6
    n_steps_kl_warmup: int | None = None
    n_epochs_kl_warmup: int | None = 400
    reduce_lr_on_plateau: bool = False
    lr_factor: float = 0.6
    lr_patience: int = 30
    lr_threshold: float = 0.0
    lr_scheduler_metric: Literal[
        "elbo_validation",
        "reconstruction_loss_validation",
        "kl_local_validation",
    ] = "elbo_validation"
    lr_min: float = 0.0
    adversarial_classifier: bool | Any = False
    scale_adversarial_loss: float | Literal["auto"] = "auto"
    compile: bool = False
    compile_kwargs: dict | None = None
    loss_kwargs: dict[str, Any] = field(default_factory=dict)

    def to_kwargs(self) -> dict[str, Any]:
        kwargs = {
            "optimizer": self.optimizer,
            "optimizer_creator": self.optimizer_creator,
            "lr": self.lr,
            "weight_decay": self.weight_decay,
            "n_steps_kl_warmup": self.n_steps_kl_warmup,
            "n_epochs_kl_warmup": self.n_epochs_kl_warmup,
            "reduce_lr_on_plateau": self.reduce_lr_on_plateau,
            "lr_factor": self.lr_factor,
            "lr_patience": self.lr_patience,
            "lr_threshold": self.lr_threshold,
            "lr_scheduler_metric": self.lr_scheduler_metric,
            "lr_min": self.lr_min,
            "adversarial_classifier": self.adversarial_classifier,
            "scale_adversarial_loss": self.scale_adversarial_loss,
            "compile": self.compile,
            "compile_kwargs": self.compile_kwargs,
        }
        kwargs.update(self.loss_kwargs)
        return kwargs


@dataclass
class SemiSupervisedTrainingPlanConfig:
    """Config for :class:`~scvi.train.SemiSupervisedTrainingPlan`."""

    classification_ratio: int = 50
    lr: float = 1e-3
    weight_decay: float = 1e-6
    n_steps_kl_warmup: int | None = None
    n_epochs_kl_warmup: int | None = 400
    reduce_lr_on_plateau: bool = False
    lr_factor: float = 0.6
    lr_patience: int = 30
    lr_threshold: float = 0.0
    lr_scheduler_metric: Literal[
        "elbo_validation",
        "reconstruction_loss_validation",
        "kl_local_validation",
    ] = "elbo_validation"
    compile: bool = False
    compile_kwargs: dict | None = None
    loss_kwargs: dict[str, Any] = field(default_factory=dict)

    def to_kwargs(self) -> dict[str, Any]:
        kwargs = {
            "classification_ratio": self.classification_ratio,
            "lr": self.lr,
            "weight_decay": self.weight_decay,
            "n_steps_kl_warmup": self.n_steps_kl_warmup,
            "n_epochs_kl_warmup": self.n_epochs_kl_warmup,
            "reduce_lr_on_plateau": self.reduce_lr_on_plateau,
            "lr_factor": self.lr_factor,
            "lr_patience": self.lr_patience,
            "lr_threshold": self.lr_threshold,
            "lr_scheduler_metric": self.lr_scheduler_metric,
            "compile": self.compile,
            "compile_kwargs": self.compile_kwargs,
        }
        kwargs.update(self.loss_kwargs)
        return kwargs


@dataclass
class SemiSupervisedAdversarialTrainingPlanConfig:
    """Config for :class:`~scvi.train.SemiSupervisedAdversarialTrainingPlan`."""

    key_adversarial: str = REGISTRY_KEYS.BATCH_KEY
    optimizer: Literal["Adam", "AdamW", "Custom"] = "Adam"
    optimizer_creator: TorchOptimizerCreator | None = None
    classification_ratio: int = 50
    lr: float = 1e-3
    weight_decay: float = 1e-6
    n_steps_kl_warmup: int | None = None
    n_epochs_kl_warmup: int | None = 400
    reduce_lr_on_plateau: bool = False
    lr_factor: float = 0.6
    lr_patience: int = 30
    lr_threshold: float = 0.0
    lr_scheduler_metric: Literal[
        "elbo_validation",
        "reconstruction_loss_validation",
        "kl_local_validation",
    ] = "elbo_validation"
    lr_min: float = 0.0
    adversarial_classifier: bool | Any = False
    scale_adversarial_loss: float | Literal["auto"] = "auto"
    loss_kwargs: dict[str, Any] = field(default_factory=dict)

    def to_kwargs(self) -> dict[str, Any]:
        kwargs = {
            "key_adversarial": self.key_adversarial,
            "optimizer": self.optimizer,
            "optimizer_creator": self.optimizer_creator,
            "classification_ratio": self.classification_ratio,
            "lr": self.lr,
            "weight_decay": self.weight_decay,
            "n_steps_kl_warmup": self.n_steps_kl_warmup,
            "n_epochs_kl_warmup": self.n_epochs_kl_warmup,
            "reduce_lr_on_plateau": self.reduce_lr_on_plateau,
            "lr_factor": self.lr_factor,
            "lr_patience": self.lr_patience,
            "lr_threshold": self.lr_threshold,
            "lr_scheduler_metric": self.lr_scheduler_metric,
            "lr_min": self.lr_min,
            "adversarial_classifier": self.adversarial_classifier,
            "scale_adversarial_loss": self.scale_adversarial_loss,
        }
        kwargs.update(self.loss_kwargs)
        return kwargs


@dataclass
class LowLevelPyroTrainingPlanConfig:
    """Config for :class:`~scvi.train.LowLevelPyroTrainingPlan`."""

    loss_fn: Any | None = None
    optim: Any | None = None
    optim_kwargs: dict | None = None
    n_steps_kl_warmup: int | None = None
    n_epochs_kl_warmup: int | None = 400
    scale_elbo: float = 1.0

    def to_kwargs(self) -> dict[str, Any]:
        return {
            "loss_fn": self.loss_fn,
            "optim": self.optim,
            "optim_kwargs": self.optim_kwargs,
            "n_steps_kl_warmup": self.n_steps_kl_warmup,
            "n_epochs_kl_warmup": self.n_epochs_kl_warmup,
            "scale_elbo": self.scale_elbo,
        }


@dataclass
class PyroTrainingPlanConfig:
    """Config for :class:`~scvi.train.PyroTrainingPlan`."""

    loss_fn: Any | None = None
    optim: Any | None = None
    optim_kwargs: dict | None = None
    n_steps_kl_warmup: int | None = None
    n_epochs_kl_warmup: int | None = 400
    scale_elbo: float = 1.0
    blocked: list | None = None

    def to_kwargs(self) -> dict[str, Any]:
        return {
            "loss_fn": self.loss_fn,
            "optim": self.optim,
            "optim_kwargs": self.optim_kwargs,
            "n_steps_kl_warmup": self.n_steps_kl_warmup,
            "n_epochs_kl_warmup": self.n_epochs_kl_warmup,
            "scale_elbo": self.scale_elbo,
            "blocked": self.blocked,
        }


@dataclass
class ClassifierTrainingPlanConfig:
    """Config for :class:`~scvi.train.ClassifierTrainingPlan`."""

    lr: float = 1e-3
    weight_decay: float = 1e-6
    eps: float = 0.01
    optimizer: Literal["Adam", "AdamW"] = "Adam"
    data_key: str = REGISTRY_KEYS.X_KEY
    labels_key: str = REGISTRY_KEYS.LABELS_KEY
    loss: Callable = torch.nn.CrossEntropyLoss

    def to_kwargs(self) -> dict[str, Any]:
        return {
            "lr": self.lr,
            "weight_decay": self.weight_decay,
            "eps": self.eps,
            "optimizer": self.optimizer,
            "data_key": self.data_key,
            "labels_key": self.labels_key,
            "loss": self.loss,
        }


@dataclass
class JaxTrainingPlanConfig:
    """Config for :class:`~scvi.train.JaxTrainingPlan`."""

    optimizer: Literal["Adam", "AdamW", "Custom"] = "Adam"
    optimizer_creator: Callable[[], Any] | None = None
    lr: float = 1e-3
    weight_decay: float = 1e-6
    eps: float = 0.01
    max_norm: float | None = None
    n_steps_kl_warmup: int | None = None
    n_epochs_kl_warmup: int | None = 400
    loss_kwargs: dict[str, Any] = field(default_factory=dict)

    def to_kwargs(self) -> dict[str, Any]:
        kwargs = {
            "optimizer": self.optimizer,
            "optimizer_creator": self.optimizer_creator,
            "lr": self.lr,
            "weight_decay": self.weight_decay,
            "eps": self.eps,
            "max_norm": self.max_norm,
            "n_steps_kl_warmup": self.n_steps_kl_warmup,
            "n_epochs_kl_warmup": self.n_epochs_kl_warmup,
        }
        kwargs.update(self.loss_kwargs)
        return kwargs


@dataclass
class TrainerConfig:
    """Config for :class:`~scvi.train.Trainer`."""

    accelerator: str | Accelerator | None = None
    devices: list[int] | str | int | None = None
    benchmark: bool = True
    check_val_every_n_epoch: int | None = None
    max_epochs: int = 400
    default_root_dir: str | None = None
    enable_checkpointing: bool = False
    checkpointing_monitor: str = "validation_loss"
    num_sanity_val_steps: int = 0
    enable_model_summary: bool = False
    early_stopping: bool = False
    early_stopping_monitor: Literal[
        "elbo_validation",
        "reconstruction_loss_validation",
        "kl_local_validation",
    ] = "elbo_validation"
    early_stopping_min_delta: float = 0.0
    early_stopping_patience: int = 45
    early_stopping_warmup_epochs: int = 0
    early_stopping_mode: Literal["min", "max"] = "min"
    enable_progress_bar: bool = True
    progress_bar_refresh_rate: int = 1
    simple_progress_bar: bool = True
    logger: Logger | None | bool = None
    log_every_n_steps: int = 10
    learning_rate_monitor: bool = False
    log_save_dir: str | None = None
    extra_kwargs: dict[str, Any] = field(default_factory=dict)

    def to_kwargs(self) -> dict[str, Any]:
        kwargs = dict(self.extra_kwargs)
        kwargs.update(
            {
                "accelerator": self.accelerator,
                "devices": self.devices,
                "benchmark": self.benchmark,
                "check_val_every_n_epoch": self.check_val_every_n_epoch,
                "max_epochs": self.max_epochs,
                "default_root_dir": self.default_root_dir,
                "enable_checkpointing": self.enable_checkpointing,
                "checkpointing_monitor": self.checkpointing_monitor,
                "num_sanity_val_steps": self.num_sanity_val_steps,
                "enable_model_summary": self.enable_model_summary,
                "early_stopping": self.early_stopping,
                "early_stopping_monitor": self.early_stopping_monitor,
                "early_stopping_min_delta": self.early_stopping_min_delta,
                "early_stopping_patience": self.early_stopping_patience,
                "early_stopping_warmup_epochs": self.early_stopping_warmup_epochs,
                "early_stopping_mode": self.early_stopping_mode,
                "enable_progress_bar": self.enable_progress_bar,
                "progress_bar_refresh_rate": self.progress_bar_refresh_rate,
                "simple_progress_bar": self.simple_progress_bar,
                "logger": self.logger,
                "log_every_n_steps": self.log_every_n_steps,
                "learning_rate_monitor": self.learning_rate_monitor,
                "log_save_dir": self.log_save_dir,
            }
        )
        return kwargs
