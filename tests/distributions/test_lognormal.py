import pytest
import torch

from scvi.distributions import Log1pNormal, LogNormal, ZeroInflatedLogNormal


def test_lognormal_distribution():
    """Test standard LogNormal distribution."""
    mu = torch.tensor([0.0, 1.0])
    sigma = torch.tensor([1.0, 0.5])

    dist = LogNormal(mu=mu, sigma=sigma, validate_args=True)

    # Test mean property (E[X] = exp(μ + σ²/2))
    expected_mean = torch.exp(mu + 0.5 * sigma**2)
    assert torch.allclose(dist.mean, expected_mean, atol=1e-6)

    # Test variance property
    sigma_sq = sigma**2
    expected_var = (torch.exp(sigma_sq) - 1) * torch.exp(2 * mu + sigma_sq)
    assert torch.allclose(dist.variance, expected_var, atol=1e-6)

    # Test sampling shapes
    s1 = dist.sample((100,))
    assert s1.shape == (100, 2)
    s2 = dist.sample(sample_shape=(4, 3))
    assert s2.shape == (4, 3, 2)

    # Test samples are positive (support is (0, ∞))
    samples = dist.sample((1000,))
    assert (samples > 0).all()

    # Test log_prob shape
    x = torch.rand(2) + 0.1  # positive values
    log_p = dist.log_prob(x)
    assert log_p.shape == (2,)

    # Test log_prob is finite for positive values
    assert torch.isfinite(log_p).all()

    # Test __repr__
    dist.__repr__()

    # Test get_normalized
    assert torch.allclose(dist.get_normalized("mu"), mu)

    # Test with scale parameter
    scale = torch.tensor([0.5, 0.5])
    dist_with_scale = LogNormal(mu=mu, sigma=sigma, scale=scale, validate_args=True)
    assert torch.allclose(dist_with_scale.get_normalized("scale"), scale)


def test_log1p_normal_distribution():
    """Test Log1pNormal distribution (log(X+1) ~ Normal)."""
    mu = torch.tensor([0.0, 1.0])
    sigma = torch.tensor([1.0, 0.5])

    dist = Log1pNormal(mu=mu, sigma=sigma, validate_args=True)

    # Test mean property (E[X] = exp(μ + σ²/2) - 1)
    expected_mean = torch.exp(mu + 0.5 * sigma**2) - 1
    assert torch.allclose(dist.mean, expected_mean, atol=1e-6)

    # Test variance property
    sigma_sq = sigma**2
    expected_var = (torch.exp(sigma_sq) - 1) * torch.exp(2 * mu + sigma_sq)
    assert torch.allclose(dist.variance, expected_var, atol=1e-6)

    # Test sampling shapes
    s1 = dist.sample((100,))
    assert s1.shape == (100, 2)
    s2 = dist.sample(sample_shape=(4, 3))
    assert s2.shape == (4, 3, 2)

    # Test samples are non-negative (support is [0, ∞))
    samples = dist.sample((1000,))
    assert (samples >= 0).all()

    # Test log_prob shape
    x = torch.rand(2)  # non-negative values
    log_p = dist.log_prob(x)
    assert log_p.shape == (2,)

    # Test log_prob at zero is finite
    log_p_zero = dist.log_prob(torch.zeros(2))
    assert torch.isfinite(log_p_zero).all()

    # Test __repr__
    dist.__repr__()

    # Test get_normalized
    assert torch.allclose(dist.get_normalized("mu"), mu)


def test_zero_inflated_lognormal_distribution():
    """Test ZeroInflatedLogNormal distribution."""
    mu = torch.tensor([0.0, 1.0])
    sigma = torch.tensor([1.0, 0.5])
    zi_logits = torch.tensor([0.0, -1.0])  # ~50% and ~27% zero inflation

    dist = ZeroInflatedLogNormal(mu=mu, sigma=sigma, zi_logits=zi_logits, validate_args=True)

    # Test zi_probs
    expected_zi_probs = torch.sigmoid(zi_logits)
    assert torch.allclose(dist.zi_probs, expected_zi_probs, atol=1e-6)

    # Test mean property (E[X] = (1 - π) * E[LogNormal])
    lognormal_mean = torch.exp(mu + 0.5 * sigma**2)
    expected_mean = (1 - expected_zi_probs) * lognormal_mean
    assert torch.allclose(dist.mean, expected_mean, atol=1e-6)

    # Test sampling shapes
    torch.manual_seed(0)
    s1 = dist.sample((100,))
    assert s1.shape == (100, 2)
    s2 = dist.sample(sample_shape=(4, 3))
    assert s2.shape == (4, 3, 2)

    # Test samples contain zeros due to zero inflation
    samples = dist.sample((1000,))
    assert (samples == 0).any()
    assert (samples >= 0).all()

    # Test log_prob shape
    x = torch.rand(2) + 0.1
    log_p = dist.log_prob(x)
    assert log_p.shape == (2,)

    # Test log_prob at zero is finite (important for ZI distributions)
    log_p_zero = dist.log_prob(torch.zeros(2))
    assert torch.isfinite(log_p_zero).all()

    # Test __repr__
    dist.__repr__()

    # Test with scale parameter
    scale = torch.tensor([0.5, 0.5])
    dist_with_scale = ZeroInflatedLogNormal(
        mu=mu, sigma=sigma, zi_logits=zi_logits, scale=scale, validate_args=True
    )
    dist_with_scale.__repr__()

    # Test different batch sizes
    size = (50, 3)
    mu_batch = torch.randn(size)
    sigma_batch = torch.rand(size) + 0.1
    zi_logits_batch = torch.randn(size)
    x_batch = torch.rand(size) + 0.1

    dist_batch = ZeroInflatedLogNormal(
        mu=mu_batch, sigma=sigma_batch, zi_logits=zi_logits_batch, validate_args=True
    )
    assert dist_batch.log_prob(x_batch).shape == size


def test_lognormal_log_prob_manual():
    """Test log_prob calculation manually for LogNormal."""
    mu = torch.tensor([0.0])
    sigma = torch.tensor([1.0])
    x = torch.tensor([1.0])

    dist = LogNormal(mu=mu, sigma=sigma)
    log_p = dist.log_prob(x)

    # Manual calculation: log p(x) = -0.5*log(2π) - log(σ) - 0.5*((log(x) - μ)/σ)² - log(x)
    y = torch.log(x)
    expected = (
        -0.5 * torch.log(2 * torch.tensor(torch.pi))
        - torch.log(sigma)
        - 0.5 * ((y - mu) / sigma) ** 2
        - torch.log(x)
    )
    assert torch.allclose(log_p, expected, atol=1e-6)


def test_log1p_normal_log_prob_manual():
    """Test log_prob calculation manually for Log1pNormal."""
    mu = torch.tensor([0.0])
    sigma = torch.tensor([1.0])
    x = torch.tensor([1.0])

    dist = Log1pNormal(mu=mu, sigma=sigma)
    log_p = dist.log_prob(x)

    # Manual calculation: log p(x) = -0.5*log(2π) - log(σ) - 0.5*((log(x+1) - μ)/σ)² - log(x+1)
    y = torch.log1p(x)
    expected = (
        -0.5 * torch.log(2 * torch.tensor(torch.pi))
        - torch.log(sigma)
        - 0.5 * ((y - mu) / sigma) ** 2
        - torch.log1p(x)
    )
    assert torch.allclose(log_p, expected, atol=1e-6)


def test_sampling_statistics():
    """Test that sample statistics match expected moments."""
    mu = torch.tensor([1.0])
    sigma = torch.tensor([0.5])

    # Test LogNormal
    dist_ln = LogNormal(mu=mu, sigma=sigma)
    samples_ln = dist_ln.sample((10000,))
    assert (samples_ln.mean(0) - dist_ln.mean).abs().mean() <= 0.5
    assert (samples_ln.var(0) - dist_ln.variance).abs().mean() <= 1.0

    # Test Log1pNormal
    dist_l1p = Log1pNormal(mu=mu, sigma=sigma)
    samples_l1p = dist_l1p.sample((10000,))
    assert (samples_l1p.mean(0) - dist_l1p.mean).abs().mean() <= 0.5

    # Test ZeroInflatedLogNormal - check zero proportion
    zi_logits = torch.tensor([0.0])  # 50% zeros expected
    dist_ziln = ZeroInflatedLogNormal(mu=mu, sigma=sigma, zi_logits=zi_logits)
    samples_ziln = dist_ziln.sample((10000,))
    zero_proportion = (samples_ziln == 0).float().mean()
    assert (zero_proportion - 0.5).abs() <= 0.05
