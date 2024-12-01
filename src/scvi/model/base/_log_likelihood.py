from __future__ import annotations

from inspect import signature
from typing import TYPE_CHECKING

import torch

if TYPE_CHECKING:
    from collections.abc import Callable, Iterator
    from typing import Any

    from torch import Tensor

    from scvi.module.base import LossOutput


def compute_elbo(
    module: Callable[[dict[str, Tensor | None], dict], tuple[Any, Any, LossOutput]],
    dataloader: Iterator[dict[str, Tensor | None]],
    return_mean: bool = True,
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
        An iterator over minibatches of data on which to compute the metric. The minibatches
        should be formatted as a dictionary of :class:`~torch.Tensor` with keys as expected by
        the ``forward`` method of ``module``.
    **kwargs
        Additional keyword arguments to pass into ``module``.
    return_mean
        If ``True``, return the mean ELBO across the dataset. If ``False``, return the ELBO for
        each cell individually.

    Returns
    -------
    The evidence lower bound (ELBO) of the data.
    """
    elbo = []
    if "full_forward_pass" in signature(module._get_inference_input).parameters:
        get_inference_input_kwargs = {"full_forward_pass": True}
    else:
        get_inference_input_kwargs = {}
    for tensors in dataloader:
        _, _, losses = module(
            tensors, **kwargs, get_inference_input_kwargs=get_inference_input_kwargs
        )
        if isinstance(losses.reconstruction_loss, dict):
            reconstruction_loss = torch.stack(list(losses.reconstruction_loss.values())).sum(dim=0)
        else:
            reconstruction_loss = losses.reconstruction_loss
        if isinstance(losses.kl_local, dict):
            kl_local = torch.stack(list(losses.kl_local.values())).sum(dim=0)
        else:
            kl_local = losses.kl_local
        elbo.append(reconstruction_loss + kl_local)

    elbo = torch.cat(elbo, dim=0)
    if return_mean:
        elbo = elbo.mean()
    return elbo


def compute_reconstruction_error(
    module: Callable[[dict[str, Tensor | None], dict], tuple[Any, Any, LossOutput]],
    dataloader: Iterator[dict[str, Tensor | None]],
    return_mean: bool = True,
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
        An iterator over minibatches of data on which to compute the metric. The minibatches
        should be formatted as a dictionary of :class:`~torch.Tensor` with keys as expected by
        the ``forward`` method of ``module``.
    return_mean
        If ``True``, return the mean reconstruction error across the dataset. If ``False``,
        return the reconstruction error for each cell individually.
    **kwargs
        Additional keyword arguments to pass into ``module``.

    Returns
    -------
    A dictionary of the reconstruction error of the data.
    """
    # Iterate once over the data and computes the reconstruction error
    if "full_forward_pass" in signature(module._get_inference_input).parameters:
        get_inference_input_kwargs = {"full_forward_pass": True}
    else:
        get_inference_input_kwargs = {}
    log_lkl = {}
    for tensors in dataloader:
        _, _, loss_output = module(
            tensors,
            loss_kwargs={"kl_weight": 1},
            get_inference_input_kwargs=get_inference_input_kwargs,
            **kwargs,
        )
        if not isinstance(loss_output.reconstruction_loss, dict):
            rec_loss_dict = {"reconstruction_loss": loss_output.reconstruction_loss}
        else:
            rec_loss_dict = loss_output.reconstruction_loss
        for key, value in rec_loss_dict.items():
            if key in log_lkl:
                log_lkl[key].append(value)
            else:
                log_lkl[key] = [value]

    for key, _ in log_lkl.items():
        log_lkl[key] = torch.cat(log_lkl[key], dim=0)
        if return_mean:
            log_lkl[key] = torch.mean(log_lkl[key])

    return log_lkl
