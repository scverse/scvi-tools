import pytest
import torch

import scvi.distributions._beta_binomial as beta_binomial_module
from scvi.distributions import BetaBinomial


def test_beta_binomial_distribution():
    alpha = torch.ones(size=(2,))
    beta = torch.ones_like(alpha)

    num_trials = 100
    total_count = torch.ones_like(alpha) * num_trials
    x = torch.randint_like(total_count, high=num_trials)

    dist_param_1 = BetaBinomial(
        total_count=total_count, alpha=alpha, beta=beta, validate_args=True
    )

    log_p_alpha_beta = dist_param_1.log_prob(x)

    mu = alpha / (alpha + beta)
    gamma = 1 / (alpha + beta + 1)

    dist_param_2 = BetaBinomial(total_count=total_count, mu=mu, gamma=gamma, validate_args=True)

    log_p_mu_gamma = dist_param_2.log_prob(x)
    assert (log_p_alpha_beta - log_p_mu_gamma).abs().max().item() <= 1e-8

    # Should fail with value outside of distribution's support
    with pytest.raises(ValueError):
        dist_param_1.log_prob(-x)

    # Should fail as no parameterization is specified
    with pytest.raises(ValueError):
        BetaBinomial(
            total_count=total_count,
        )

    # Should fail as two full parameterizations are provided
    with pytest.raises(ValueError):
        BetaBinomial(total_count=total_count, alpha=alpha, beta=beta, mu=mu, gamma=gamma)

    # Should fail without a complete parameterization 1 or 2
    with pytest.raises(ValueError):
        BetaBinomial(total_count=total_count, alpha=alpha, gamma=gamma)


def test_sample_cpu_detour_forced(monkeypatch):
    """Forces the CPU detour the way a torch without MPS '_sample_dirichlet'/'binomial' kernels
    would."""
    monkeypatch.setattr(beta_binomial_module, "_needs_cpu_detour", lambda on_mps, op: True)

    alpha = torch.ones(64, 8) * 2
    beta = torch.ones(64, 8) * 2
    total_count = torch.ones(64, 8) * 100
    dist = BetaBinomial(total_count=total_count, alpha=alpha, beta=beta)

    samples = dist.sample((16,))
    assert samples.shape == (16, 64, 8)
    assert (samples >= 0).all()
    assert (samples <= 100).all()
    assert samples.device == alpha.device


@pytest.mark.skipif(not torch.backends.mps.is_available(), reason="requires an MPS device")
def test_sample_stays_on_mps():
    alpha = torch.ones(64, 8, device="mps") * 2
    beta = torch.ones(64, 8, device="mps") * 2
    total_count = torch.ones(64, 8, device="mps") * 100
    dist = BetaBinomial(total_count=total_count, alpha=alpha, beta=beta)

    samples = dist.sample((16,))
    assert samples.device.type == "mps"
    assert samples.shape == (16, 64, 8)
    assert (samples >= 0).all()
    assert (samples <= 100).all()
