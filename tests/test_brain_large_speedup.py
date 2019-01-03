#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""Test improved calculation of log_zinb_forward for better performance on brain_large."""

from scipy.sparse import random
import torch
import torch.nn.functional as F

batch_size = 128
num_genes = 720

x = torch.from_numpy(random(batch_size, num_genes,
                            density=.5).A).type(torch.float32)
mu = torch.rand(batch_size, num_genes)
theta = torch.rand(batch_size, num_genes)
pi = torch.rand(batch_size, num_genes)
eps = 1e-3


def test_log_zinb_positive():
    """Test that the new method to compute log_zinb_positive is the same as the existing."""

    def existing_method(x, mu, theta, pi, eps=1e-8):
        case_zero = (F.softplus((- pi + theta * torch.log(theta + eps) - theta * torch.log(theta + mu + eps))) -
                     F.softplus(-pi))

        case_non_zero = - pi - F.softplus(-pi) + \
            theta * torch.log(theta + eps) - \
            theta * torch.log(theta + mu + eps) + \
            x * torch.log(mu + eps) - \
            x * torch.log(theta + mu + eps) + \
            torch.lgamma(x + theta) - \
            torch.lgamma(theta) - \
            torch.lgamma(x + 1)

        res = torch.mul((x < eps).type(torch.float32), case_zero) + \
            torch.mul((x > eps).type(torch.float32), case_non_zero)

        return torch.sum(res, dim=-1)

    def new_method(x, mu, theta, pi, eps=1e-8):
        softplus_pi = F.softplus(-pi)
        log_theta_eps = torch.log(theta + eps)
        log_theta_mu_eps = torch.log(theta + mu + eps)
        pi_theta_log = - pi + theta * (log_theta_eps - log_theta_mu_eps)

        case_zero = F.softplus(pi_theta_log) - softplus_pi
        mul_case_zero = torch.mul((x < eps).type(torch.float32), case_zero)

        case_non_zero = - softplus_pi + \
            pi_theta_log + \
            x * (torch.log(mu + eps) - log_theta_mu_eps) + \
            torch.lgamma(x + theta) - \
            torch.lgamma(theta) - \
            torch.lgamma(x + 1)

        res = mul_case_zero + torch.mul((x > eps).type(torch.float32), case_non_zero)

        return torch.sum(res, dim=-1)

    existing_likelihood = existing_method(x, mu, theta, pi)
    new_likelihood = new_method(x, mu, theta, pi)

    diff = torch.abs(existing_likelihood - new_likelihood)

    assert (diff < eps).all()


def test_log_nb_forward():
    """Test that the new method to compute log_nb_positive is the same as the existing."""

    def existing_method(x, mu, theta, eps=1e-8):
        res = theta * torch.log(theta + eps) - theta * torch.log(theta + mu + eps) + x * torch.log(
            mu + eps) - x * torch.log(theta + mu + eps) + torch.lgamma(x + theta) - torch.lgamma(
            theta) - torch.lgamma(x + 1)

        return torch.sum(res, dim=-1)

    def new_method(x, mu, theta, eps=1e-8):
        res = theta * (torch.log(theta + eps) - torch.log(theta + mu + eps)) + \
            x * (torch.log(mu + eps) - torch.log(theta + mu + eps)) + \
            torch.lgamma(x + theta) - \
            torch.lgamma(theta) - \
            torch.lgamma(x + 1)

        return torch.sum(res, dim=-1)

    existing_likelihood = existing_method(x, mu, theta, pi)
    new_likelihood = new_method(x, mu, theta, pi)

    diff = torch.abs(existing_likelihood - new_likelihood)

    assert (diff < eps).all()
