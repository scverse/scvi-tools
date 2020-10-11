"""File for computing log likelihood of the data"""
import numpy as np
import torch

from torch import logsumexp
from torch.distributions import Beta, Normal

from scvi import _CONSTANTS


def compute_elbo(vae, data_loader, **kwargs):
    """
    Computes the ELBO.

    The ELBO is the reconstruction error + the KL divergences
    between the variational distributions and the priors.
    It differs from the marginal log likelihood.
    Specifically, it is a lower bound on the marginal log likelihood
    plus a term that is constant with respect to the variational distribution.
    It still gives good insights on the modeling of the data, and is fast to compute.
    """
    # Iterate once over the data and compute the elbo
    elbo = 0
    for i_batch, tensors in enumerate(data_loader):
        # sample_batch = tensors[_CONSTANTS.X_KEY]
        # local_l_mean = tensors[_CONSTANTS.LOCAL_L_MEAN_KEY]
        # local_l_var = tensors[_CONSTANTS.LOCAL_L_VAR_KEY]
        # batch_index = tensors[_CONSTANTS.BATCH_KEY]
        # labels = tensors[_CONSTANTS.LABELS_KEY]

        # kl_divergence_global (scalar) should be common across all batches after training
        # reconst_loss, kl_divergence, kl_divergence_global = vae(
        #     sample_batch,
        #     local_l_mean,
        #     local_l_var,
        #     batch_index=batch_index,
        #     y=labels,
        #     **kwargs
        # )
        # elbo += torch.sum(reconst_loss + kl_divergence).item()
        _, losses = vae(tensors)

        reconstruction_losses = losses["reconstruction_losses"]
        kls_local = losses["kl_local"]

        if isinstance(reconstruction_losses, dict):
            reconstruction_loss = 0.0
            for value in reconstruction_losses.values():
                reconstruction_loss += value
        else:
            reconstruction_loss = reconstruction_losses

        if isinstance(kls_local, dict):
            kl_local = 0.0
            for kl in kls_local.values():
                kl_local += kl
        else:
            kl_local = kls_local

        elbo += torch.sum(reconstruction_loss + kl_local).item()

    kl_globals = losses["kl_global"]
    if isinstance(kl_globals, dict):
        kl_global = 0.0
        for kl in kl_globals.values():
            kl_global += kl
    else:
        kl_global = kl_globals

    n_samples = len(data_loader.indices)
    # elbo += kl_divergence_global
    return elbo / n_samples


def compute_reconstruction_error(vae, data_loader, **kwargs):
    """
    Computes log p(x/z), which is the reconstruction error.

    Differs from the marginal log likelihood, but still gives good
    insights on the modeling of the data, and is fast to compute.
    """
    # Iterate once over the data and computes the reconstruction error
    log_lkl = 0
    for i_batch, tensors in enumerate(data_loader):
        # sample_batch = tensors[_CONSTANTS.X_KEY]
        # batch_index = tensors[_CONSTANTS.BATCH_KEY]
        # labels = tensors[_CONSTANTS.LABELS_KEY]

        # Distribution parameters
        # outputs = vae.inference(sample_batch, batch_index, labels, **kwargs)
        # px_r = outputs["px_r"]
        # px_rate = outputs["px_rate"]
        # px_dropout = outputs["px_dropout"]
        # bernoulli_params = outputs.get("bernoulli_params", None)

        # should call forward and take the first term
        # outputs = vae.inference(sample_batch, batch_index, labels, **kwargs)
        loss_kwargs = dict(kl_weight=1)
        _, losses = vae(tensors, loss_kwargs=loss_kwargs)
        reconstruction_loss = losses["reconstruction_loss"]

        # this if for TotalVI
        if isinstance(reconstruction_loss, dict):
            recon_loss = 0
            for value in recon_loss.values():
                recon_loss += value
            reconstruction_loss = recon_loss

        # Reconstruction loss
        # reconst_loss = vae.get_reconstruction_loss(
        #     sample_batch,
        #     px_rate,
        #     px_r,
        #     px_dropout,
        #     bernoulli_params=bernoulli_params,
        #     **kwargs
        # )

        # log_lkl += torch.sum(reconst_loss).item()
        log_lkl += torch.sum(reconstruction_loss).item()

    n_samples = len(data_loader.indices)
    return log_lkl / n_samples


def compute_marginal_log_likelihood_scvi(vae, data_loader, n_samples_mc=100):
    """
    Computes a biased estimator for log p(x), which is the marginal log likelihood.

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
    for i_batch, tensors in enumerate(data_loader):
        sample_batch = tensors[_CONSTANTS.X_KEY]
        local_l_mean = tensors[_CONSTANTS.LOCAL_L_MEAN_KEY]
        local_l_var = tensors[_CONSTANTS.LOCAL_L_VAR_KEY]
        batch_index = tensors[_CONSTANTS.BATCH_KEY]
        labels = tensors[_CONSTANTS.LABELS_KEY]

        to_sum = torch.zeros(sample_batch.size()[0], n_samples_mc)

        for i in range(n_samples_mc):

            # Distribution parameters and sampled variables
            outputs = vae.inference(sample_batch, batch_index, labels)
            # px_r = outputs["px_r"]
            # px_rate = outputs["px_rate"]
            # px_dropout = outputs["px_dropout"]
            qz_m = outputs["qz_m"]
            qz_v = outputs["qz_v"]
            z = outputs["z"]
            ql_m = outputs["ql_m"]
            ql_v = outputs["ql_v"]
            library = outputs["library"]

            # Reconstruction Loss
            # reconst_loss = vae.get_reconstruction_loss(
            #     sample_batch, px_rate, px_r, px_dropout
            # )
            reconst_loss = vae.get_reconstruction_loss(sample_batch, outputs)

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

    n_samples = len(data_loader.indices)
    # The minus sign is there because we actually look at the negative log likelihood
    return -log_lkl / n_samples


def compute_marginal_log_likelihood_autozi(autozivae, data_loader, n_samples_mc=100):
    """
    Computes a biased estimator for log p(x), which is the marginal log likelihood.

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

        for i_batch, tensors in enumerate(data_loader):
            sample_batch = tensors[_CONSTANTS.X_KEY]
            local_l_mean = tensors[_CONSTANTS.LOCAL_L_MEAN_KEY]
            local_l_var = tensors[_CONSTANTS.LOCAL_L_VAR_KEY]
            batch_index = tensors[_CONSTANTS.BATCH_KEY]
            labels = tensors[_CONSTANTS.LABELS_KEY]

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
    n_samples = len(data_loader.indices)
    # The minus sign is there because we actually look at the negative log likelihood
    return -log_lkl / n_samples
