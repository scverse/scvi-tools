import numpy as np
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


@pytest.mark.parametrize("sparse_format", [None, "csr_matrix"])
def test_scrna_raw_counts_properties_matches_dense_reference(sparse_format):
    from scvi import REGISTRY_KEYS
    from scvi.data import synthetic_iid
    from scvi.model import SCVI
    from scvi.model._utils import scrna_raw_counts_properties

    adata = synthetic_iid(sparse_format=sparse_format)
    SCVI.setup_anndata(adata, batch_key="batch")
    manager = SCVI(adata).adata_manager

    idx1 = np.array([0, 3, 7, 10, 25])
    idx2 = np.array([1, 2, 4, 8, 16, 32])
    props = scrna_raw_counts_properties(manager, idx1, idx2)

    x = manager.get_from_registry(REGISTRY_KEYS.X_KEY)
    x = np.asarray(x.todense()) if sparse_format is not None else np.asarray(x)
    scaling = 1e4 / x.sum(axis=1, keepdims=True)
    for i, idx in ((1, idx1), (2, idx2)):
        np.testing.assert_allclose(props[f"raw_mean{i}"], x[idx].mean(axis=0), rtol=1e-6)
        np.testing.assert_allclose(
            props[f"non_zeros_proportion{i}"], (x[idx] != 0).mean(axis=0), rtol=1e-6
        )
        np.testing.assert_allclose(
            props[f"raw_normalized_mean{i}"], (x[idx] * scaling[idx]).mean(axis=0), rtol=1e-6
        )
