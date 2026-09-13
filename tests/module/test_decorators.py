import pytest
import torch
from torch import nn

from scvi.module.base import _decorators as decorators_module
from scvi.module.base import auto_move_data


class _Model(nn.Module):
    def __init__(self):
        super().__init__()
        self.linear = nn.Linear(4, 2)

    @auto_move_data
    def forward(self, x: torch.Tensor, extra: dict | None = None):
        return self.linear(x), extra


class _ParameterlessModel(nn.Module):
    @auto_move_data
    def forward(self, x: torch.Tensor):
        return x


def test_auto_move_data_targets_the_parameter_device(monkeypatch):
    requested = []
    original = decorators_module._move_data_to_device

    def spy(batch, device):
        requested.append(device)
        return original(batch, device)

    monkeypatch.setattr(decorators_module, "_move_data_to_device", spy)

    model = _Model()
    x = torch.zeros(3, 4)
    extra = {"a": torch.zeros(1), "b": "not a tensor"}

    model.train()
    out, extra_out = model(x, extra=extra)
    assert out.shape == (3, 2)
    assert extra_out is extra
    # in training mode nothing is moved
    assert requested == []

    model.eval()
    out, extra_out = model(x, extra=extra)
    assert out.shape == (3, 2)
    assert extra_out["b"] == "not a tensor"
    # args and kwargs are each moved once, to the device the parameters are on
    assert requested == [model.linear.weight.device] * 2

    if torch.cuda.is_available():
        model.cuda()
        out, extra_out = model(x, extra=extra)
        assert out.device.type == "cuda"
        assert extra_out["a"].device.type == "cuda"


def test_auto_move_data_without_parameters():
    model = _ParameterlessModel().eval()
    with pytest.raises(RuntimeError, match="no parameters"):
        model(torch.zeros(2))
