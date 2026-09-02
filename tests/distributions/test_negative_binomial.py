import pytest
import torch

from scvi.distributions import (
    NegativeBinomial,
    NegativeBinomialMixture,
    ZeroInflatedNegativeBinomial,
)
from scvi.distributions import _negative_binomial as nb_module
from scvi.distributions import _utils as dist_utils
from scvi.distributions._negative_binomial import (
    _gamma,
    log_nb_positive,
    log_zinb_positive,
)
from scvi.distributions._utils import _mps_supports, _needs_cpu_detour


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
    with pytest.warns(UserWarning, match="The value argument must be within the support"):
        dist1.log_prob(-x)  # ensures neg values raise warning
    with pytest.warns(UserWarning, match="The value argument must be within the support"):
        dist2.log_prob(0.5 * x)  # ensures float values raise warning

    # test with no scale
    dist1 = ZeroInflatedNegativeBinomial(mu=mu, theta=theta, zi_logits=pi, validate_args=True)
    dist2 = NegativeBinomial(mu=mu, theta=theta, validate_args=True)
    dist1.__repr__()
    dist2.__repr__()
    assert dist1.log_prob(x).shape == size
    assert dist2.log_prob(x).shape == size

    with pytest.warns(UserWarning, match="The value argument must be within the support"):
        dist1.log_prob(-x)
    with pytest.warns(UserWarning, match="The value argument must be within the support"):
        dist2.log_prob(0.5 * x)


@pytest.mark.skipif(not torch.backends.mps.is_available(), reason="requires an MPS device")
def test_sampling_stays_on_mps():
    theta = 100.0 + torch.rand(size=(64, 8), device="mps")
    mu = 15.0 * torch.ones_like(theta)
    dists = [
        NegativeBinomial(mu=mu, theta=theta),
        ZeroInflatedNegativeBinomial(mu=mu, theta=theta, zi_logits=torch.randn_like(theta)),
        NegativeBinomialMixture(
            mu1=mu, mu2=2.0 * mu, theta1=theta, mixture_logits=torch.zeros_like(theta)
        ),
    ]
    for dist in dists:
        samples = dist.sample((16,))
        assert samples.device.type == "mps"
        assert samples.shape == (16, 64, 8)
        assert (samples >= 0).all()


@pytest.mark.skipif(not torch.backends.mps.is_available(), reason="requires an MPS device")
def test_gamma_stays_on_mps_when_the_kernel_exists():
    # The CPU detour is conditional on the running torch build, not hardcoded,
    # so assert against the same probe rather than against a version.
    theta = 100.0 + torch.rand(size=(4,), device="mps")
    mu = 15.0 * torch.ones_like(theta)
    expected = "mps" if _mps_supports(torch._standard_gamma) else "cpu"
    assert _gamma(theta, mu, on_mps=True).concentration.device.type == expected


@pytest.mark.parametrize(
    ("on_mps", "kernel_present", "expected"),
    [(False, False, False), (False, True, False), (True, True, False), (True, False, True)],
)
def test_cpu_detour_is_taken_only_for_mps_without_a_kernel(
    monkeypatch, on_mps, kernel_present, expected
):
    monkeypatch.setattr(dist_utils, "_mps_supports", lambda op: kernel_present)
    assert _needs_cpu_detour(on_mps, torch.poisson) is expected


def test_mps_support_probe_answers_without_a_device():
    # The probe builds its own mps tensor inside the try, so a machine with no MPS
    # gets False from the same code path rather than an exception.
    supported = _mps_supports(torch.poisson)
    assert isinstance(supported, bool)
    if not torch.backends.mps.is_available():
        assert supported is False


def test_sampling_returns_the_device_it_was_given():
    theta = 100.0 + torch.rand(size=(8, 3))
    mu = 15.0 * torch.ones_like(theta)
    for dist in (
        NegativeBinomial(mu=mu, theta=theta),
        NegativeBinomialMixture(
            mu1=mu, mu2=2.0 * mu, theta1=theta, mixture_logits=torch.zeros_like(theta)
        ),
    ):
        assert dist.sample((4,)).device == mu.device


def test_cpu_detour_path_still_samples_correctly(monkeypatch):
    # Force the fallback the way a torch without the kernels would, so the detour
    # itself is exercised on any machine rather than only on Apple silicon.
    monkeypatch.setattr(nb_module, "_needs_cpu_detour", lambda on_mps, op: True)
    theta = 100.0 + torch.rand(size=(2,))
    mu = 15.0 * torch.ones_like(theta)
    samples = NegativeBinomial(mu=mu, theta=theta).sample((1000,))
    assert samples.shape == (1000, 2)
    assert samples.device == mu.device
    assert (samples.mean(0) - mu).abs().mean() <= 1e0
    assert _gamma(theta, mu, on_mps=True).concentration.device.type == "cpu"
