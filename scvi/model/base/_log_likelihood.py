"""File for computing log likelihood of the data."""

import numpy as np
import torch


def compute_elbo(vae, data_loader, feed_labels=True, return_mean=True, **kwargs):
    """Computes the ELBO.

    The ELBO is the reconstruction error + the KL divergences
    between the variational distributions and the priors.
    It differs from the marginal log likelihood.
    Specifically, it is a lower bound on the marginal log likelihood
    plus a term that is constant with respect to the variational distribution.
    It still gives good insights on the modeling of the data, and is fast to compute.
    """
    # Iterate once over the data and compute the elbo
    if return_mean:
        elbo = 0
    else:
        elbo = np.array([])
    for tensors in data_loader:
        _, _, scvi_loss = vae(tensors, **kwargs)

        recon_loss = np.sum([np.array(i) for i in scvi_loss.reconstruction_loss.values()], axis=0)
        kl_local = np.sum([np.array(i) for i in scvi_loss.kl_local.values()], axis=0)

        if return_mean:
            elbo += (recon_loss + kl_local).sum(0).item()
        else:
            elbo = np.concatenate((elbo, recon_loss + kl_local), axis=0)

    kl_global = np.sum([np.array(i) for i in scvi_loss.kl_global.values()], axis=0)
    n_samples = len(data_loader.indices)
    if return_mean:
        elbo += kl_global
        return elbo / n_samples
    else:
        return elbo + kl_global / n_samples


# do each one
def compute_reconstruction_error(vae, data_loader, return_mean=True, **kwargs):
    """Computes log p(x/z), which is the reconstruction error.

    Differs from the marginal log likelihood, but still gives good
    insights on the modeling of the data, and is fast to compute.
    """
    # Iterate once over the data and computes the reconstruction error
    log_lkl = {}
    for tensors in data_loader:
        loss_kwargs = {"kl_weight": 1}
        _, _, losses = vae(tensors, loss_kwargs=loss_kwargs)
        if not isinstance(losses.reconstruction_loss, dict):
            rec_loss_dict = {"reconstruction_loss": losses.reconstruction_loss}
        else:
            rec_loss_dict = losses.reconstruction_loss
        for key, value in rec_loss_dict.items():
            if return_mean:
                if key in log_lkl:
                    if return_mean:
                        log_lkl[key] += torch.sum(value).item()
                    else:
                        log_lkl[key].append(value)
                else:
                    if return_mean:
                        log_lkl[key] = torch.sum(value).item()
                    else:
                        log_lkl[key] = value

    n_samples = len(data_loader.indices)
    for key, _ in log_lkl.items():
        if return_mean:
            log_lkl[key] = log_lkl[key] / n_samples
        else:
            log_lkl[key] = torch.cat(log_lkl[key], dim=0)

        log_lkl[key] = -log_lkl[key]
    return log_lkl
