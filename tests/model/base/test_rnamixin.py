import pytest
import torch

import scvi.model.base._rnamixin as rnamixin_module
from scvi.data import synthetic_iid
from scvi.model import SCVI


def _make_trained_model(gene_likelihood: str = "nb"):
    adata = synthetic_iid()
    SCVI.setup_anndata(adata, batch_key="batch")
    model = SCVI(adata, gene_likelihood=gene_likelihood)
    model.train(1, check_val_every_n_epoch=1, train_size=0.5)
    return model


def test_posterior_predictive_sample_gamma_cpu_detour_forced(monkeypatch):
    """Forces the CPU detour the way a torch without an MPS '_standard_gamma' kernel would."""
    monkeypatch.setattr(rnamixin_module, "_needs_cpu_detour", lambda on_mps, op: True)

    model = _make_trained_model(gene_likelihood="nb")
    sample = model.posterior_predictive_sample(n_samples=2)
    assert sample.shape == (*model.adata.shape, 2)


@pytest.mark.skipif(not torch.backends.mps.is_available(), reason="requires an MPS device")
def test_posterior_predictive_sample_on_mps():
    adata = synthetic_iid()
    SCVI.setup_anndata(adata, batch_key="batch")
    model = SCVI(adata, gene_likelihood="nb")
    model.train(1, check_val_every_n_epoch=1, train_size=0.5, accelerator="mps")
    sample = model.posterior_predictive_sample(n_samples=2)
    assert sample.shape == (*adata.shape, 2)
