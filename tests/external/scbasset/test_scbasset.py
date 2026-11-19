import numpy as np
import pytest
import torch
import torchmetrics

from scvi import REGISTRY_KEYS
from scvi.data import synthetic_iid
from scvi.distributions._utils import _mps_supports
from scvi.external import SCBASSET

_DNA_CODE_KEY = "code"


def _get_adata(sparse_format: str | None = None):
    dataset1 = synthetic_iid(batch_size=100, sparse_format=sparse_format).transpose()
    dataset1.X = (dataset1.X > 0).astype(float)
    dataset1.obsm[_DNA_CODE_KEY] = np.random.randint(0, 3, size=(dataset1.n_obs, 1344))
    return dataset1


@pytest.mark.optional
def test_scbasset():
    adata = _get_adata()
    SCBASSET.setup_anndata(
        adata,
        dna_code_key=_DNA_CODE_KEY,
    )
    model = SCBASSET(adata)
    model.train(max_epochs=2, early_stopping=True)
    model.get_latent_representation()


@pytest.mark.optional
def test_scbasset_batch():
    adata = _get_adata()
    SCBASSET.setup_anndata(
        adata,
        dna_code_key=_DNA_CODE_KEY,
        batch_key="batch",
    )
    model = SCBASSET(adata)
    model.train(max_epochs=1)
    model.get_latent_representation()
    assert hasattr(model.module, "batch_ids")


def test_scbasset_auroc_matches_cpu_on_large_batch():
    # torchmetrics' binary AUROC is wrong on MPS once the input reaches 65,536 elements
    # (upstream F.pad bug), so the loss must route it through the CPU there. Elsewhere
    # this runs on CPU and checks the same value, so it holds on every platform. The
    # upstream failure is non-deterministic (it can return the right value by chance),
    # so the device check is what reliably catches a regression on MPS.
    device = "mps" if _mps_supports(torch.sigmoid) else "cpu"
    adata = synthetic_iid(batch_size=300, n_genes=20).transpose()  # 600 cells
    adata.X = (adata.X > 0).astype(np.float32)
    adata.obsm[_DNA_CODE_KEY] = np.random.randint(0, 3, size=(adata.n_obs, 1344))
    SCBASSET.setup_anndata(adata, dna_code_key=_DNA_CODE_KEY)
    module = SCBASSET(adata).module.to(device)

    generator = torch.Generator().manual_seed(0)
    logits = torch.randn(128, adata.n_obs, generator=generator)  # 76,800 > 65,536
    target = (torch.rand(128, adata.n_obs, generator=generator) < 0.05).float()
    expected = torchmetrics.functional.auroc(
        torch.sigmoid(logits).ravel(), target.int().ravel(), task="binary"
    )

    loss_output = module.loss(
        {REGISTRY_KEYS.X_KEY: target.to(device)},
        {},
        {"reconstruction_logits": logits.to(device)},
    )
    auroc = loss_output.extra_metrics["auroc"]
    assert auroc.item() == pytest.approx(expected.item(), abs=1e-6)
    assert auroc.device.type == ("cpu" if device == "mps" else device)
