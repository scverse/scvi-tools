"""File for computing log likelihood of the data"""

import torch
import torch.nn.functional as F
from torch.distributions import Normal
import numpy as np


def compute_log_likelihood(vae, data_loader):
    # Iterate once over the data_loader and computes the total log_likelihood
    log_lkl = 0
    for i_batch, tensors in enumerate(data_loader):
        sample_batch, local_l_mean, local_l_var, batch_index, labels = tensors
        reconst_loss, kl_divergence = vae(sample_batch, local_l_mean, local_l_var, batch_index=batch_index,
                                          y=labels)
        log_lkl += torch.sum(reconst_loss).item()
    n_samples = (len(data_loader.dataset)
                 if not (hasattr(data_loader, 'sampler') and hasattr(data_loader.sampler, 'indices')) else
                 len(data_loader.sampler.indices))
    return log_lkl / n_samples


def compute_tighter_log_likelihood(vae, data_loader, n_samples_mc=100):
    # Uses MC sampling to compute a tighter lower bound
    log_lkl = 0
    for i_batch, tensors in enumerate(data_loader):
        sample_batch, local_l_mean, local_l_var, batch_index, labels = tensors
        batch_log_lkl = compute_tighter_log_likelihood_sample(vae, sample_batch, local_l_mean,
                                                              local_l_var, batch_index,
                                                              labels, n_samples_mc=n_samples_mc)
        log_lkl += torch.sum(batch_log_lkl).item()
    n_samples = (len(data_loader.dataset)
                 if not (hasattr(data_loader, 'sampler') and hasattr(data_loader.sampler, 'indices')) else
                 len(data_loader.sampler.indices))
    return log_lkl / n_samples


def logsumexp(inputs, dim=None, keepdim=False):
    """Numerically stable logsumexp.

    Args:
        inputs: A Variable with any shape.
        dim: An integer.
        keepdim: A boolean.

    Returns:
        Equivalent of log(sum(exp(inputs), dim=dim, keepdim=keepdim)).
    """
    # For a 1-D array x (any array along a single dimension),
    # log sum exp(x) = s + log sum exp(x - s)
    # with s = max(x) being a common choice.
    if dim is None:
        inputs = inputs.view(-1)
        dim = 0
    s, _ = torch.max(inputs, dim=dim, keepdim=True)
    outputs = s + (inputs - s).exp().sum(dim=dim, keepdim=True).log()
    if not keepdim:
        outputs = outputs.squeeze(dim)
    return outputs


def compute_tighter_log_likelihood_sample(vae, sample_batch, local_l_mean, local_l_var, batch_index,
                                          labels, n_samples_mc=100):
    # Uses MC sampling to compute a tighter lower bound
    x = torch.log(1 + sample_batch)
    to_sum = torch.zeros(sample_batch.size()[0], n_samples_mc)
    for i in range(n_samples_mc):
        qz_m, qz_v, z = vae.z_encoder(x, labels)
        reconst_loss, kl_divergence = vae(sample_batch, local_l_mean,
                                          local_l_var,
                                          batch_index=batch_index,
                                          y=labels)
        p_z = Normal(torch.zeros_like(z), torch.ones_like(z)).log_prob(z).sum(dim=-1)
        p_x_z = reconst_loss
        q_z_x = Normal(qz_m, qz_v).log_prob(z).sum(dim=-1)
        to_sum[:, i] = p_z + p_x_z - q_z_x - np.log(n_samples_mc)

    return logsumexp(to_sum, dim=-1)


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

    case_zero = (F.softplus((- pi + theta * torch.log(theta + eps) - theta * torch.log(theta + mu + eps)))
                 - F.softplus(-pi))

    case_non_zero = - pi - F.softplus(-pi) + theta * torch.log(theta + eps) - theta * torch.log(
        theta + mu + eps) + x * torch.log(mu + eps) - x * torch.log(theta + mu + eps) + torch.lgamma(
        x + theta + eps) - torch.lgamma(theta + eps) - torch.lgamma(x + 1)

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
