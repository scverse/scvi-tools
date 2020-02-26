"""File for computing log likelihood of the data"""

import numpy as np
import torch
import torch.nn.functional as F
from torch import logsumexp
from torch.distributions import Normal, Beta


def compute_elbo(vae, posterior, **kwargs):
    """ Computes the ELBO.

    The ELBO is the reconstruction error + the KL divergences
    between the variational distributions and the priors.
    It differs from the marginal log likelihood.
    Specifically, it is a lower bound on the marginal log likelihood
    plus a term that is constant with respect to the variational distribution.
    It still gives good insights on the modeling of the data, and is fast to compute.
    """
    # Iterate once over the posterior and compute the elbo
    elbo = 0
    for i_batch, tensors in enumerate(posterior):
        sample_batch, local_l_mean, local_l_var, batch_index, labels = tensors[
            :5
        ]  # general fish case
        # kl_divergence_global (scalar) should be common across all batches after training
        reconst_loss, kl_divergence, kl_divergence_global = vae(
            sample_batch,
            local_l_mean,
            local_l_var,
            batch_index=batch_index,
            y=labels,
            **kwargs
        )
        elbo += torch.sum(reconst_loss + kl_divergence).item()
    n_samples = len(posterior.indices)
    elbo += kl_divergence_global
    return elbo / n_samples


def compute_reconstruction_error(vae, posterior, **kwargs):
    """ Computes log p(x/z), which is the reconstruction error.

    Differs from the marginal log likelihood, but still gives good
    insights on the modeling of the data, and is fast to compute.
    """
    # Iterate once over the posterior and computes the reconstruction error
    log_lkl = 0
    for i_batch, tensors in enumerate(posterior):
        sample_batch, local_l_mean, local_l_var, batch_index, labels = tensors[
            :5
        ]  # general fish case

        # Distribution parameters
        outputs = vae.inference(sample_batch, batch_index, labels, **kwargs)
        px_r = outputs["px_r"]
        px_rate = outputs["px_rate"]
        px_dropout = outputs["px_dropout"]
        bernoulli_params = outputs.get("bernoulli_params", None)

        # Reconstruction loss
        reconst_loss = vae.get_reconstruction_loss(
            sample_batch,
            px_rate,
            px_r,
            px_dropout,
            bernoulli_params=bernoulli_params,
            **kwargs
        )

        log_lkl += torch.sum(reconst_loss).item()
    n_samples = len(posterior.indices)
    return log_lkl / n_samples


def compute_marginal_log_likelihood_scvi(vae, posterior, n_samples_mc=100):
    """ Computes a biased estimator for log p(x), which is the marginal log likelihood.

    Despite its bias, the estimator still converges to the real value
    of log p(x) when n_samples_mc (for Monte Carlo) goes to infinity
    (a fairly high value like 100 should be enough)
    Due to the Monte Carlo sampling, this method is not as computationally efficient
    as computing only the reconstruction loss
    """
    if vae.latent_distribution == "ln":
        raise NotImplementedError

    # Uses MC sampling to compute a tighter lower bound on log p(x)
    log_lkl = 0
    for i_batch, tensors in enumerate(posterior):
        sample_batch, local_l_mean, local_l_var, batch_index, labels = tensors
        to_sum = torch.zeros(sample_batch.size()[0], n_samples_mc)

        for i in range(n_samples_mc):

            # Distribution parameters and sampled variables
            outputs = vae.inference(sample_batch, batch_index, labels)
            px_r = outputs["px_r"]
            px_rate = outputs["px_rate"]
            px_dropout = outputs["px_dropout"]
            qz_m = outputs["qz_m"]
            qz_v = outputs["qz_v"]
            z = outputs["z"]
            ql_m = outputs["ql_m"]
            ql_v = outputs["ql_v"]
            library = outputs["library"]

            # Reconstruction Loss
            reconst_loss = vae.get_reconstruction_loss(
                sample_batch, px_rate, px_r, px_dropout
            )

            # Log-probabilities
            p_l = Normal(local_l_mean, local_l_var.sqrt()).log_prob(library).sum(dim=-1)
            p_z = (
                Normal(torch.zeros_like(qz_m), torch.ones_like(qz_v))
                .log_prob(z)
                .sum(dim=-1)
            )
            p_x_zl = -reconst_loss
            q_z_x = Normal(qz_m, qz_v.sqrt()).log_prob(z).sum(dim=-1)
            q_l_x = Normal(ql_m, ql_v.sqrt()).log_prob(library).sum(dim=-1)

            to_sum[:, i] = p_z + p_l + p_x_zl - q_z_x - q_l_x

        batch_log_lkl = logsumexp(to_sum, dim=-1) - np.log(n_samples_mc)
        log_lkl += torch.sum(batch_log_lkl).item()

    n_samples = len(posterior.indices)
    # The minus sign is there because we actually look at the negative log likelihood
    return -log_lkl / n_samples


def compute_marginal_log_likelihood_autozi(autozivae, posterior, n_samples_mc=100):
    """ Computes a biased estimator for log p(x), which is the marginal log likelihood.

    Despite its bias, the estimator still converges to the real value
    of log p(x) when n_samples_mc (for Monte Carlo) goes to infinity
    (a fairly high value like 100 should be enough)
    Due to the Monte Carlo sampling, this method is not as computationally efficient
    as computing only the reconstruction loss
    """
    # Uses MC sampling to compute a tighter lower bound on log p(x)
    log_lkl = 0
    to_sum = torch.zeros((n_samples_mc,))
    alphas_betas = autozivae.get_alphas_betas(as_numpy=False)
    alpha_prior = alphas_betas["alpha_prior"]
    alpha_posterior = alphas_betas["alpha_posterior"]
    beta_prior = alphas_betas["beta_prior"]
    beta_posterior = alphas_betas["beta_posterior"]

    for i in range(n_samples_mc):

        bernoulli_params = autozivae.sample_from_beta_distribution(
            alpha_posterior, beta_posterior
        )

        for i_batch, tensors in enumerate(posterior):
            sample_batch, local_l_mean, local_l_var, batch_index, labels = tensors

            # Distribution parameters and sampled variables
            outputs = autozivae.inference(sample_batch, batch_index, labels)
            px_r = outputs["px_r"]
            px_rate = outputs["px_rate"]
            px_dropout = outputs["px_dropout"]
            qz_m = outputs["qz_m"]
            qz_v = outputs["qz_v"]
            z = outputs["z"]
            ql_m = outputs["ql_m"]
            ql_v = outputs["ql_v"]
            library = outputs["library"]

            # Reconstruction Loss
            bernoulli_params_batch = autozivae.reshape_bernoulli(
                bernoulli_params, batch_index, labels
            )
            reconst_loss = autozivae.get_reconstruction_loss(
                sample_batch, px_rate, px_r, px_dropout, bernoulli_params_batch
            )

            # Log-probabilities
            p_l = Normal(local_l_mean, local_l_var.sqrt()).log_prob(library).sum(dim=-1)
            p_z = (
                Normal(torch.zeros_like(qz_m), torch.ones_like(qz_v))
                .log_prob(z)
                .sum(dim=-1)
            )
            p_x_zld = -reconst_loss
            q_z_x = Normal(qz_m, qz_v.sqrt()).log_prob(z).sum(dim=-1)
            q_l_x = Normal(ql_m, ql_v.sqrt()).log_prob(library).sum(dim=-1)

            batch_log_lkl = torch.sum(p_x_zld + p_l + p_z - q_z_x - q_l_x, dim=0)
            to_sum[i] += batch_log_lkl

        p_d = Beta(alpha_prior, beta_prior).log_prob(bernoulli_params).sum()
        q_d = Beta(alpha_posterior, beta_posterior).log_prob(bernoulli_params).sum()

        to_sum[i] += p_d - q_d

    log_lkl = logsumexp(to_sum, dim=-1).item() - np.log(n_samples_mc)
    n_samples = len(posterior.indices)
    # The minus sign is there because we actually look at the negative log likelihood
    return -log_lkl / n_samples


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
        theta = theta.view(
            1, theta.size(0)
        )  # In this case, we reshape theta for broadcasting

    softplus_pi = F.softplus(-pi)  #  uses log(sigmoid(x)) = -softplus(-x)
    log_theta_eps = torch.log(theta + eps)
    log_theta_mu_eps = torch.log(theta + mu + eps)
    pi_theta_log = -pi + theta * (log_theta_eps - log_theta_mu_eps)

    case_zero = F.softplus(pi_theta_log) - softplus_pi
    mul_case_zero = torch.mul((x < eps).type(torch.float32), case_zero)

    case_non_zero = (
        -softplus_pi
        + pi_theta_log
        + x * (torch.log(mu + eps) - log_theta_mu_eps)
        + torch.lgamma(x + theta)
        - torch.lgamma(theta)
        - torch.lgamma(x + 1)
    )
    mul_case_non_zero = torch.mul((x > eps).type(torch.float32), case_non_zero)

    res = mul_case_zero + mul_case_non_zero

    return res


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
        theta = theta.view(
            1, theta.size(0)
        )  # In this case, we reshape theta for broadcasting

    log_theta_mu_eps = torch.log(theta + mu + eps)

    res = (
        theta * (torch.log(theta + eps) - log_theta_mu_eps)
        + x * (torch.log(mu + eps) - log_theta_mu_eps)
        + torch.lgamma(x + theta)
        - torch.lgamma(theta)
        - torch.lgamma(x + 1)
    )

    return res


def log_mixture_nb(x, mu_1, mu_2, theta_1, theta_2, pi, eps=1e-8):
    """
    Note: All inputs should be torch Tensors
    log likelihood (scalar) of a minibatch according to a mixture nb model.
    pi is the probability to be in the first component.

    For totalVI, the first component should be background.

    Variables:
    mu1: mean of the first negative binomial component (has to be positive support) (shape: minibatch x genes)
    theta1: first inverse dispersion parameter (has to be positive support) (shape: minibatch x genes)
    mu2: mean of the second negative binomial (has to be positive support) (shape: minibatch x genes)
    theta2: second inverse dispersion parameter (has to be positive support) (shape: minibatch x genes)
        If None, assume one shared inverse dispersion parameter.
    eps: numerical stability constant
    """
    if theta_2 is not None:
        log_nb_1 = log_nb_positive(x, mu_1, theta_1)
        log_nb_2 = log_nb_positive(x, mu_2, theta_2)
    # this is intended to reduce repeated computations
    else:
        theta = theta_1
        if theta.ndimension() == 1:
            theta = theta.view(
                1, theta.size(0)
            )  # In this case, we reshape theta for broadcasting

        log_theta_mu_1_eps = torch.log(theta + mu_1 + eps)
        log_theta_mu_2_eps = torch.log(theta + mu_2 + eps)
        lgamma_x_theta = torch.lgamma(x + theta)
        lgamma_theta = torch.lgamma(theta)
        lgamma_x_plus_1 = torch.lgamma(x + 1)

        log_nb_1 = (
            theta * (torch.log(theta + eps) - log_theta_mu_1_eps)
            + x * (torch.log(mu_1 + eps) - log_theta_mu_1_eps)
            + lgamma_x_theta
            - lgamma_theta
            - lgamma_x_plus_1
        )
        log_nb_2 = (
            theta * (torch.log(theta + eps) - log_theta_mu_2_eps)
            + x * (torch.log(mu_2 + eps) - log_theta_mu_2_eps)
            + lgamma_x_theta
            - lgamma_theta
            - lgamma_x_plus_1
        )

    logsumexp = torch.logsumexp(torch.stack((log_nb_1, log_nb_2 - pi)), dim=0)
    softplus_pi = F.softplus(-pi)

    log_mixture_nb = logsumexp - softplus_pi

    return log_mixture_nb
