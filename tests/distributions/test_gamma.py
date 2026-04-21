import torch

from scvi.distributions import ZeroInflatedGamma


def test_zero_inflated_gamma_distribution():
    """Test ZeroInflatedGamma distribution."""
    concentration = torch.tensor([2.0, 3.0])
    rate = torch.tensor([1.0, 2.0])
    zi_logits = torch.tensor([0.0, -1.0])  # ~50% and ~27% zero inflation

    dist = ZeroInflatedGamma(
        concentration=concentration, rate=rate, zi_logits=zi_logits, validate_args=True
    )

    # Test zi_probs
    expected_zi_probs = torch.sigmoid(zi_logits)
    assert torch.allclose(dist.zi_probs, expected_zi_probs, atol=1e-6)

    # Test mean property (E[X] = (1 - π) * E[Gamma])
    gamma_mean = concentration / rate
    expected_mean = (1 - expected_zi_probs) * gamma_mean
    assert torch.allclose(dist.mean, expected_mean, atol=1e-6)

    # Test variance property
    gamma_var = concentration / (rate**2)
    expected_var = (1 - expected_zi_probs) * (gamma_var + expected_zi_probs * gamma_mean**2)
    assert torch.allclose(dist.variance, expected_var, atol=1e-6)

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
    dist_with_scale = ZeroInflatedGamma(
        concentration=concentration,
        rate=rate,
        zi_logits=zi_logits,
        scale=scale,
        validate_args=True,
    )
    dist_with_scale.__repr__()
    assert torch.allclose(dist_with_scale.get_normalized("scale"), scale)
    assert torch.allclose(dist_with_scale.get_normalized("concentration"), concentration)
    assert torch.allclose(dist_with_scale.get_normalized("rate"), rate)

    # Test different batch sizes
    size = (50, 3)
    concentration_batch = torch.rand(size) + 0.5
    rate_batch = torch.rand(size) + 0.5
    zi_logits_batch = torch.randn(size)
    x_batch = torch.rand(size) + 0.1

    dist_batch = ZeroInflatedGamma(
        concentration=concentration_batch,
        rate=rate_batch,
        zi_logits=zi_logits_batch,
        validate_args=True,
    )
    assert dist_batch.log_prob(x_batch).shape == size


def test_sampling_statistics_zi_gamma():
    """Test that sample statistics match expected moments for ZeroInflatedGamma."""
    concentration = torch.tensor([5.0])
    rate = torch.tensor([2.0])

    # Test ZeroInflatedGamma - check zero proportion
    zi_logits = torch.tensor([0.0])  # 50% zeros expected
    dist_zig = ZeroInflatedGamma(concentration=concentration, rate=rate, zi_logits=zi_logits)
    samples_zig = dist_zig.sample((10000,))
    zero_proportion = (samples_zig == 0).float().mean()
    assert (zero_proportion - 0.5).abs() <= 0.05


def test_zi_gamma_high_zero_inflation():
    """Test ZeroInflatedGamma with high zero inflation rate."""
    concentration = torch.tensor([2.0])
    rate = torch.tensor([1.0])
    zi_logits = torch.tensor([2.0])  # ~88% zero inflation

    dist = ZeroInflatedGamma(concentration=concentration, rate=rate, zi_logits=zi_logits)
    samples = dist.sample((1000,))

    # Check that most samples are zeros
    zero_proportion = (samples == 0).float().mean()
    expected_zi_prob = torch.sigmoid(zi_logits).item()
    assert (zero_proportion - expected_zi_prob).abs() <= 0.05


def test_zi_gamma_low_zero_inflation():
    """Test ZeroInflatedGamma with low zero inflation rate."""
    concentration = torch.tensor([2.0])
    rate = torch.tensor([1.0])
    zi_logits = torch.tensor([-3.0])  # ~5% zero inflation

    dist = ZeroInflatedGamma(concentration=concentration, rate=rate, zi_logits=zi_logits)
    samples = dist.sample((1000,))

    # Check that few samples are zeros
    zero_proportion = (samples == 0).float().mean()
    expected_zi_prob = torch.sigmoid(zi_logits).item()
    assert (zero_proportion - expected_zi_prob).abs() <= 0.05
