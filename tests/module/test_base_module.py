import pytest
from torch import nn

from scvi.module.base._base_module import BaseModuleClass


def test_device_targets_the_parameter_device():
    module = BaseModuleClass()
    module.add_module("linear", nn.Linear(4, 2))
    assert module.device == module.linear.weight.device


def test_device_without_parameters_raises():
    module = BaseModuleClass()
    with pytest.raises(RuntimeError, match="no parameters"):
        _ = module.device
