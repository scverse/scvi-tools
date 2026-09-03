import pytest
import torch

from scvi.model._utils import parse_device_args


@pytest.mark.skipif(not torch.backends.mps.is_available(), reason="requires an MPS device")
def test_auto_accelerator_falls_back_to_cpu_when_mps_available():
    with pytest.warns(UserWarning, match="automatically set to `cpu`"):
        _, _, device = parse_device_args(accelerator="auto", devices="auto", return_device="torch")
    assert device.type == "cpu"


@pytest.mark.skipif(not torch.backends.mps.is_available(), reason="requires an MPS device")
def test_explicit_mps_accelerator_warns_about_backend_caveats():
    with pytest.warns(UserWarning, match="Results will not be bit-identical"):
        _, _, device = parse_device_args(accelerator="mps", devices="auto", return_device="torch")
    assert device.type == "mps"
