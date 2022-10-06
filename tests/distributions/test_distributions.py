import numpy as np
import pytest
import torch

from scvi.distributions import NegativeBinomial, ZeroInflatedNegativeBinomial
from scvi.distributions._negative_binomial import log_nb_positive, log_zinb_positive
from scvi.model._metrics import unsupervised_clustering_accuracy

use_gpu = True


def test_deprecated_munkres():
    y = np.array([0, 1, 0, 1, 0, 1, 1, 1])
    y_pred = np.array([0, 0, 0, 0, 1, 1, 1, 1])
    reward, assignment = unsupervised_clustering_accuracy(y, y_pred)
    assert reward == 0.625
    assert (assignment == np.array([[0, 0], [1, 1]])).all()

    y = np.array([1, 1, 2, 2, 0, 0, 3, 3])
    y_pred = np.array([1, 1, 2, 2, 3, 3, 0, 0])
    reward, assignment = unsupervised_clustering_accuracy(y, y_pred)
    assert reward == 1.0
    assert (assignment == np.array([[0, 3], [1, 1], [2, 2], [3, 0]])).all()


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
    pi = torch.randn_like(theta)
    x = torch.randint_like(mu, high=20)
    dist1 = ZeroInflatedNegativeBinomial(
        mu=mu, theta=theta, zi_logits=pi, validate_args=True
    )
    dist2 = NegativeBinomial(mu=mu, theta=theta, validate_args=True)
    assert dist1.log_prob(x).shape == size
    assert dist2.log_prob(x).shape == size

    with pytest.raises(ValueError):
        ZeroInflatedNegativeBinomial(
            mu=-mu, theta=theta, zi_logits=pi, validate_args=True
        )
    with pytest.warns(UserWarning):
        dist1.log_prob(-x)  # ensures neg values raise warning
    with pytest.warns(UserWarning):
        dist2.log_prob(0.5 * x)  # ensures float values raise warning
