"""File for computing log likelihood of the data"""

import torch

from functions.gamma import Lgamma


def log_zinb_positive(x, mu, theta, pi, eps=1e-8):
    """
    Note: All inputs are torch Tensors
    log likelihood (scalar) of a minibatch according to a zinb model.
    Notes:
    We parametrize the bernoulli using the logits, hence the softplus functions appearing

    Variables:
    mu: mean of the negative binomial (has to be positive support) (shape: minibatch x genes)
    theta: inverse dispersion parameter (has to be positive support) (shape: minibatch x genes)
    pi: logit of the dropout parameter (real support) (shape: minibatch x genes)
    eps: numerical stability constant
    """
    def softplus(x):
        return torch.log(1 + torch.exp(x))

    case_zero = softplus((- pi + theta * torch.log(theta + eps) - theta * torch.log(theta + mu + eps)))
    - softplus(-pi)

    case_non_zero = - pi - softplus(-pi) + theta * torch.log(theta + eps) - theta * torch.log(
        theta + mu + eps) + x * torch.log(mu + eps) - x * torch.log(theta + mu + eps) + Lgamma()(
        x + theta) - Lgamma()(theta.resize(1, theta.size(0))) - Lgamma()(x + 1)

    mask = x.clone()
    mask[mask < eps] = 1
    mask[mask >= eps] = 0
    res = torch.mul(mask, case_zero) + torch.mul(1 - mask, case_non_zero)
    return torch.sum(res, dim=-1)


def log_nb_positive(x, mu, theta, eps=1e-8):
    """
    Note: All inputs should be torch Tensors
    log likelihood (scalar) of a minibatch according to a nb model.

    Variables:
    mu: mean of the negative binomial (has to be positive support) (shape: minibatch x genes)
    theta: inverse dispersion parameter (has to be positive support) (shape: minibatch x genes)
    eps: numerical stability constant
    """
    res = theta * torch.log(theta + eps) - theta * torch.log(theta + mu + eps) + x * torch.log(
        mu + eps) - x * torch.log(theta + mu + eps) + Lgamma()(x + theta) - Lgamma()(
        theta.resize(1, theta.size(0))) - Lgamma()(
        x + 1)
    return torch.sum(res)
