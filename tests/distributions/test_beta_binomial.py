import pytest
import torch

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
