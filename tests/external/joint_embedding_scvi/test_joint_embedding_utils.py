import torch

from scvi.external.joint_embedding_scvi._utils import (
    binomial_split,
    cross_correlation_loss,
    sample_thinning_probs,
    variance_loss,
)


def test_binomial_split_conserves_counts():
    x = torch.randint(0, 50, (32, 10)).float()
    x1, x2 = binomial_split(x, p=0.5)
    assert torch.all(x1 >= 0)
    assert torch.all(x2 >= 0)
    assert torch.allclose(x1 + x2, x)


def test_binomial_split_per_cell_probs():
    x = torch.randint(0, 50, (8, 5)).float()
    p = torch.rand(8, 1)
    x1, x2 = binomial_split(x, p=p)
    assert torch.allclose(x1 + x2, x)


def test_binomial_split_rounds_non_integer_counts():
    # Non-integer counts must not raise; they are rounded before thinning.
    x = torch.full((4, 3), 3.7)
    x1, x2 = binomial_split(x, p=0.5)
    assert torch.allclose(x1 + x2, torch.full((4, 3), 4.0))


def test_sample_thinning_probs_range_and_shape():
    x = torch.randint(0, 100, (16, 20)).float()
    p = sample_thinning_probs(x, min_library_size=10.0)
    assert p.shape == (16, 1)
    assert torch.all(p >= 0.0)
    assert torch.all(p <= 1.0)


def test_cross_correlation_loss_identity_components():
    z = torch.randn(64, 8)
    comp = cross_correlation_loss(z, z, lambda_off_diag=0.01, return_components=True)
    # z vs itself -> near-perfect invariance (diagonal ~1)
    assert comp["invariance"].item() < 1e-3
    # total decomposes exactly into invariance + lambda * redundancy
    assert torch.allclose(comp["total"], comp["invariance"] + 0.01 * comp["redundancy"])


def test_variance_loss_collapsed_vs_spread():
    collapsed = torch.zeros(32, 5)
    spread = torch.randn(2000, 5) * 2.0
    assert variance_loss(collapsed, target_std=1.0).item() > 0.9
    assert variance_loss(spread, target_std=1.0).item() < 0.1


def test_losses_finite_on_single_sample_batch():
    # Single-row minibatches must stay finite (biased std), not NaN.
    z = torch.randn(1, 5)
    assert torch.isfinite(variance_loss(z))
    assert torch.isfinite(cross_correlation_loss(z, z))
