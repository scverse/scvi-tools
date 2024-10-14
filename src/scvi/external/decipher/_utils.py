import numpy as np
import pyro.poutine as poutine
import torch


def predictive_log_likelihood(decipher_module, batch, n_samples=5):
    """
    Calculate the predictive log-likelihood for a Decipher module.

    This function performs multiple runs through the dataloader to obtain
    an empirical estimate of the predictive log-likelihood. It calculates the
    log-likelihood for each run and returns the average. The beta parameter
    of the Decipher module is temporarily modified and restored even if an
    exception occurs.

    Parameters
    ----------
    decipher_module : PyroBaseModuleClass
        The Decipher module to evaluate.
    batch : torch.Tensor
        Batch of data to compute the log-likelihood for.
    n_samples : int, optional
        Number of passes through the dataloader (default is 5).

    Returns
    -------
    float
        The average estimated predictive log-likelihood across multiple runs.
    """
    log_weights = []
    old_beta = decipher_module.beta
    decipher_module.beta = 1.0
    try:
        for _ in range(n_samples):
            guide_trace = poutine.trace(decipher_module.guide).get_trace(batch)
            model_trace = poutine.trace(
                poutine.replay(decipher_module.model, trace=guide_trace)
            ).get_trace(batch)
            log_weights.append(model_trace.log_prob_sum() - guide_trace.log_prob_sum())

    finally:
        decipher_module.beta = old_beta

    log_z = torch.logsumexp(torch.tensor(log_weights) - np.log(n_samples), 0)
    return log_z.item()
