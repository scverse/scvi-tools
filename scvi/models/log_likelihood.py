"""File for computing log likelihood of the data"""

import numpy as np
import torch
import torch.nn.functional as F
from torch import logsumexp
from torch.distributions import Normal


def compute_log_likelihood(vae, posterior, **kwargs):
    """ Computes log p(x/z), which is the reconstruction error .
        Differs from the marginal log likelihood, but still gives good
        insights on the modeling of the data, and is fast to compute
    """
    # Iterate once over the posterior and computes the total log_likelihood
    log_lkl = 0
    for i_batch, tensors in enumerate(posterior):
        sample_batch, local_l_mean, local_l_var, batch_index, labels = tensors[:5]  # general fish case
        reconst_loss, kl_divergence = vae(sample_batch, local_l_mean, local_l_var, batch_index=batch_index,
                                          y=labels, **kwargs)
        log_lkl += torch.sum(reconst_loss).item()
    n_samples = len(posterior.indices)
    return log_lkl / n_samples


def compute_marginal_log_likelihood(vae, posterior, n_samples_mc=100):
    """ Computes a biased estimator for log p(x), which is the marginal log likelihood.
        Despite its bias, the estimator still converges to the real value
        of log p(x) when n_samples_mc (for Monte Carlo) goes to infinity
        (a fairly high value like 100 should be enough)
        Due to the Monte Carlo sampling, this method is not as computationally efficient
        as computing only the reconstruction loss
    """
    # Uses MC sampling to compute a tighter lower bound
    log_lkl = 0
    for i_batch, tensors in enumerate(posterior):
        sample_batch, local_l_mean, local_l_var, batch_index, labels = tensors
        x = torch.log(1 + sample_batch)
        to_sum = torch.zeros(sample_batch.size()[0], n_samples_mc)
        for i in range(n_samples_mc):
            qz_m, qz_v, z = vae.z_encoder(x, labels)
            reconst_loss, kl_divergence = vae(sample_batch, local_l_mean,
                                              local_l_var,
                                              batch_index=batch_index,
                                              y=labels)
            p_z = Normal(torch.zeros_like(qz_m), torch.ones_like(qz_v)).log_prob(z).sum(dim=-1)
            p_x_z = - reconst_loss
            q_z_x = Normal(qz_m, qz_v.sqrt()).log_prob(z).sum(dim=-1)
            to_sum[:, i] = p_z + p_x_z - q_z_x
        batch_log_lkl = logsumexp(to_sum, dim=-1) - np.log(n_samples_mc)
        log_lkl += torch.sum(batch_log_lkl).item()
    n_samples = len(posterior.indices)
    # The minus sign is there because we actually look at the negative log likelihood
    return - log_lkl / n_samples


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

    # theta is the dispersion rate. If .ndimension() == 1, it is shared for all cells (regardless of batch or labels)
    if theta.ndimension() == 1:
        theta = theta.view(1, theta.size(0))  # In this case, we reshape theta for broadcasting

    case_zero = (F.softplus((- pi + theta * torch.log(theta + eps) - theta * torch.log(theta + mu + eps))) -
                 F.softplus(-pi))

    case_non_zero = - pi - F.softplus(-pi) + theta * torch.log(theta + eps) - theta * torch.log(
        theta + mu + eps) + x * torch.log(mu + eps) - x * torch.log(theta + mu + eps) + torch.lgamma(
        x + theta) - torch.lgamma(theta) - torch.lgamma(x + 1)

    res = torch.mul((x < eps).type(torch.float32), case_zero) + torch.mul((x > eps).type(torch.float32), case_non_zero)
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
    if theta.ndimension() == 1:
        theta = theta.view(1, theta.size(0))  # In this case, we reshape theta for broadcasting

    res = theta * torch.log(theta + eps) - theta * torch.log(theta + mu + eps) + x * torch.log(
        mu + eps) - x * torch.log(theta + mu + eps) + torch.lgamma(x + theta) - torch.lgamma(
        theta) - torch.lgamma(x + 1)
    return torch.sum(res, dim=-1)
