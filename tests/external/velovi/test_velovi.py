import numpy as np
import pytest
import torch

import scvi.external.velovi._module as velovi_module_module
from scvi.data import synthetic_iid
from scvi.external.velovi import VELOVI
from scvi.external.velovi._module import VELOVAE


@pytest.mark.optional
def test_velovi():
    n_latent = 5
    adata = synthetic_iid()
    adata.layers["spliced"] = adata.X.copy()
    adata.layers["unspliced"] = adata.X.copy()
    VELOVI.setup_anndata(adata, unspliced_layer="unspliced", spliced_layer="spliced")
    model = VELOVI(adata, n_latent=n_latent)
    model.train(1, check_val_every_n_epoch=1, train_size=0.5)
    model.get_latent_representation()
    model.get_velocity()
    model.get_latent_time()
    model.get_state_assignment()
    model.get_expression_fit()
    model.get_directional_uncertainty()
    model.get_permutation_scores(labels_key="labels")

    _ = model.history

    # tests __repr__
    print(model)


def _make_velovae(n_input: int) -> VELOVAE:
    return VELOVAE(
        n_input=n_input,
        switch_spliced=np.ones(n_input, dtype=np.float32),
        switch_unspliced=np.ones(n_input, dtype=np.float32),
    )


def test_generative_dirichlet_cpu_detour_forced(monkeypatch):
    """Forces the CPU detour the way a torch without an MPS '_sample_dirichlet' kernel would."""
    monkeypatch.setattr(velovi_module_module, "_needs_cpu_detour", lambda on_mps, op: True)

    n_input = 50
    module = _make_velovae(n_input)
    z = torch.randn(8, module.n_latent)
    gamma, beta, alpha, alpha_1, lambda_alpha = module._get_rates()

    outputs = module.generative(z, gamma, beta, alpha, alpha_1, lambda_alpha)
    assert outputs["px_pi"].shape == (8, n_input, 4)
    assert torch.allclose(outputs["px_pi"].sum(dim=-1), torch.ones(8, n_input), atol=1e-5)


@pytest.mark.skipif(not torch.backends.mps.is_available(), reason="requires an MPS device")
def test_generative_dirichlet_stays_on_mps():
    n_input = 50
    module = _make_velovae(n_input).to("mps")
    z = torch.randn(8, module.n_latent, device="mps")
    gamma, beta, alpha, alpha_1, lambda_alpha = module._get_rates()

    outputs = module.generative(z, gamma, beta, alpha, alpha_1, lambda_alpha)
    assert outputs["px_pi"].device.type == "mps"
    assert outputs["px_pi"].shape == (8, n_input, 4)
