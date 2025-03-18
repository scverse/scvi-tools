import pytest
import torch

from scvi.distributions import NegativeBinomial, ZeroInflatedNegativeBinomial
from scvi.distributions._negative_binomial import log_nb_positive, log_zinb_positive


def test_zinb_distribution():
    theta = 100.0 + torch.rand(size=(2,))
    mu = 15.0 * torch.ones_like(theta)
    pi = torch.randn_like(theta)
    x = torch.randint_like(mu, high=20)
    log_p_ref = log_zinb_positive(x, mu, theta, pi)

    dist = ZeroInflatedNegativeBinomial(mu=mu, theta=theta, zi_logits=pi)
    log_p_zinb = dist.log_prob(x)
    assert (log_p_ref - log_p_zinb).abs().max().item() <= 1e-8

    torch.manual_seed(0)
    s1 = dist.sample((100,))
    assert s1.shape == (100, 2)
    s2 = dist.sample(sample_shape=(4, 3))
    assert s2.shape == (4, 3, 2)

    log_p_ref = log_nb_positive(x, mu, theta)
    dist = NegativeBinomial(mu=mu, theta=theta)
    log_p_nb = dist.log_prob(x)
    assert (log_p_ref - log_p_nb).abs().max().item() <= 1e-8

    s1 = dist.sample((1000,))
    assert s1.shape == (1000, 2)
    assert (s1.mean(0) - mu).abs().mean() <= 1e0
    assert (s1.std(0) - (mu + mu * mu / theta) ** 0.5).abs().mean() <= 1e0

    size = (50, 3)
    theta = 100.0 + torch.rand(size=size)
    mu = 15.0 * torch.ones_like(theta)
    scales = mu / mu.sum(-1, keepdim=True)
    pi = torch.randn_like(theta)
    x = torch.randint_like(mu, high=20)
    dist1 = ZeroInflatedNegativeBinomial(
        mu=mu, theta=theta, zi_logits=pi, scale=scales, validate_args=True
    )
    dist2 = NegativeBinomial(mu=mu, theta=theta, scale=scales, validate_args=True)
    assert dist1.log_prob(x).shape == size
    assert dist2.log_prob(x).shape == size

    with pytest.raises(ValueError):
        ZeroInflatedNegativeBinomial(
            mu=-mu, theta=theta, zi_logits=pi, scale=scales, validate_args=True
        )
    with pytest.warns(UserWarning):
        dist1.log_prob(-x)  # ensures neg values raise warning
    with pytest.warns(UserWarning):
        dist2.log_prob(0.5 * x)  # ensures float values raise warning

    # test with no scale
    dist1 = ZeroInflatedNegativeBinomial(mu=mu, theta=theta, zi_logits=pi, validate_args=True)
    dist2 = NegativeBinomial(mu=mu, theta=theta, validate_args=True)
    dist1.__repr__()
    dist2.__repr__()
    assert dist1.log_prob(x).shape == size
    assert dist2.log_prob(x).shape == size

    with pytest.warns(UserWarning):
        dist1.log_prob(-x)
    with pytest.warns(UserWarning):
        dist2.log_prob(0.5 * x)
