import pytest
import torch

from scvi.model._utils import parse_device_args


@pytest.mark.skipif(not torch.backends.mps.is_available(), reason="requires an MPS device")
def test_auto_accelerator_falls_back_to_cpu_when_mps_available(monkeypatch):
    # Regardless of whether this suite is itself running under SCVI_ALLOW_MPS_AUTO (e.g. under
    # test_mps.yaml), this test is specifically about the default (override-unset) behavior.
    monkeypatch.delenv("SCVI_ALLOW_MPS_AUTO", raising=False)
    with pytest.warns(UserWarning, match="automatically set to `cpu`"):
        _, _, device = parse_device_args(accelerator="auto", devices="auto", return_device="torch")
    assert device.type == "cpu"


@pytest.mark.skipif(not torch.backends.mps.is_available(), reason="requires an MPS device")
def test_explicit_mps_accelerator_warns_about_backend_caveats():
    with pytest.warns(UserWarning, match="Results will not be bit-identical"):
        _, _, device = parse_device_args(accelerator="mps", devices="auto", return_device="torch")
    assert device.type == "mps"


@pytest.mark.skipif(not torch.backends.mps.is_available(), reason="requires an MPS device")
def test_auto_accelerator_resolves_to_mps_when_env_override_set(monkeypatch):
    """SCVI_ALLOW_MPS_AUTO opts a whole test run into `auto` behaving like it does on CUDA
    machines, without changing the library's real default. Used by test_mps.yaml so the
    dedicated MPS CI job actually exercises `mps` for calls that don't pass an accelerator,
    instead of only the handful of tests that consume the --accelerator/--devices fixtures."""
    monkeypatch.setenv("SCVI_ALLOW_MPS_AUTO", "1")
    with pytest.warns(UserWarning, match="Results will not be bit-identical"):
        _, _, device = parse_device_args(accelerator="auto", devices="auto", return_device="torch")
    assert device.type == "mps"
