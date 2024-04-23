from __future__ import annotations

from collections.abc import Iterator
from typing import Any, Callable

from torch import Tensor

from scvi.module.base import LossOutput


def compute_elbo(
    module: Callable[[dict[str, Tensor | None], dict], tuple[Any, Any, LossOutput]],
    dataloader: Iterator[dict[str, Tensor | None]],
    **kwargs,
) -> float:
    """Compute the evidence lower bound (ELBO) on the data.

    The ELBO is the reconstruction error plus the Kullback-Leibler (KL) divergences between the
    variational distributions and the priors. It is different from the marginal log-likelihood;
    specifically, it is a lower bound on the marginal log-likelihood plus a term that is constant
    with respect to the variational distribution. It still gives good insights on the modeling of
    the data and is fast to compute.

    Parameters
    ----------
    module
        A callable (can be a :class:`~torch.nn.Module` instance) that takes a dictionary of
        :class:`~torch.Tensor`s as input and returns a tuple of three elements, where the last
        element is an instance of :class:`~scvi.module.base.LossOutput`.
    dataloader
        An iterator that returns a dictionary of :class:`~torch.Tensor` instances for each
        minibatch formatted as expected by the ``forward`` method of ``vae``.
    **kwargs
        Additional keyword arguments to pass into ``module``.

    Returns
    -------
    The evidence lower bound (ELBO) of the data.
    """
    elbo = 0.0
    for tensors in dataloader:
        _, _, loss_output = module(tensors, **kwargs)
        elbo += (loss_output.reconstruction_loss_sum + loss_output.kl_local_sum).item()

    kl_global = loss_output.kl_global_sum
    n_samples = len(dataloader.indices)
    elbo += kl_global

    return elbo / n_samples


def compute_reconstruction_error(
    module: Callable[[dict[str, Tensor | None], dict], tuple[Any, Any, LossOutput]],
    dataloader: Iterator[dict[str, Tensor | None]],
    **kwargs,
) -> dict[str, float]:
    """Compute the reconstruction error on the data.

    The reconstruction error is the negative log likelihood of the data given the latent
    variables. It is different from the marginal log-likelihood, but still gives good insights on
    the modeling of the data and is fast to compute.

    Parameters
    ----------
    module
        A callable (can be a :class:`~torch.nn.Module` instance) that takes a dictionary of
        :class:`~torch.Tensor`s as input and returns a tuple of three elements, where the last
        element is an instance of :class:`~scvi.module.base.LossOutput`.
    dataloader
        An iterator that returns a dictionary of :class:`~torch.Tensor` instances for each
        minibatch formatted as expected by the ``forward`` method of ``vae``.
    **kwargs
        Additional keyword arguments to pass into ``module``.

    Returns
    -------
    A dictionary of the reconstruction error of the data.
    """
    log_likelihoods = {}
    for tensors in dataloader:
        _, _, loss_output = module(tensors, loss_kwargs={"kl_weight": 1}, **kwargs)
        rec_losses: dict[str, Tensor] | Tensor = loss_output.reconstruction_loss
        if not isinstance(rec_losses, dict):
            rec_losses = {"reconstruction_loss": rec_losses}

        for key, value in rec_losses.items():
            log_likelihoods[key] = log_likelihoods.get(key, 0.0) + value.sum().item()

    n_samples = len(dataloader.indices)

    return {key: -(value / n_samples) for key, value in log_likelihoods.items()}
