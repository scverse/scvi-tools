import pytest

from scvi.autotune import TunerManager
from scvi.model import SCVI


def test_tuner_manager_init():
    manager = TunerManager(SCVI)
    assert hasattr(manager, "_model_cls")
    assert hasattr(manager, "_defaults")
    assert hasattr(manager, "_registry")

    registry = manager._registry
    assert "tunables" in registry
    assert "metrics" in registry


def test_tuner_manager_basic_validation():
    manager = TunerManager(SCVI)

    # invalid params should raise an exception
    with pytest.raises(ValueError):
        manager._validate_search_space({"not_a_param": None}, False)

    # search space does not change with `use_defaults == False
    search_space = manager._validate_search_space({"n_hidden": None}, False)
    assert search_space == {"n_hidden": None}

    # invalid metrics should raise an exception
    with pytest.raises(ValueError):
        manager._validate_metrics("not_a_metric", [])
