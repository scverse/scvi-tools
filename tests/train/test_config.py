from __future__ import annotations

from types import MappingProxyType

import pytest

from scvi.train._config import TrainerConfig, TrainingPlanConfig, merge_kwargs


def test_merge_kwargs_prefers_overrides():
    config = TrainingPlanConfig(lr=1e-3, weight_decay=1e-2, loss_kwargs={"foo": "bar"})
    overrides = {"lr": 1e-4, "baz": 7}

    merged = merge_kwargs(config, overrides, name="plan")

    assert merged["lr"] == 1e-4
    assert merged["weight_decay"] == 1e-2
    assert merged["foo"] == "bar"
    assert merged["baz"] == 7


def test_merge_kwargs_accepts_mapping():
    mapping = MappingProxyType({"lr": 0.2})

    merged = merge_kwargs(mapping, None, name="plan")

    assert merged["lr"] == 0.2


def test_trainer_config_merge():
    config = TrainerConfig(check_val_every_n_epoch=2, extra_kwargs={"strategy": "ddp"})
    merged = merge_kwargs(config, {"check_val_every_n_epoch": 5}, name="trainer")

    assert merged["check_val_every_n_epoch"] == 5
    assert merged["strategy"] == "ddp"


def test_merge_kwargs_rejects_bad_config():
    class BadConfig:
        def to_kwargs(self):
            return ["not", "a", "dict"]

    with pytest.raises(TypeError):
        merge_kwargs(BadConfig(), None, name="plan")
