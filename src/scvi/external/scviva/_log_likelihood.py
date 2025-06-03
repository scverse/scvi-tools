from __future__ import annotations

from typing import TYPE_CHECKING

import torch

if TYPE_CHECKING:
    from collections.abc import Callable, Iterator
    from typing import Any

    from torch import Tensor

    from ._module import NicheLossOutput


def compute_composition_error(
    module: Callable[[dict[str, Tensor | None], dict], tuple[Any, Any, NicheLossOutput]],
    dataloader: Iterator[dict[str, Tensor | None]],
    return_mean: bool = True,
    **kwargs,
) -> float | Tensor[float]:
    """Compute the composition prediction error on the data.

    The  error is the negative log likelihood of the data (alpha) given the latent
    variables.

    Parameters
    ----------
    module
        A callable (can be a :class:`~torch.nn.Module` instance) that takes a dictionary of
        :class:`~torch.Tensor`s as input and returns a tuple of three elements, where the last
        element is an instance of :class:`~scvi.external.scviva._module.NicheLossOutput`.
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
    The composition prediction error on the data.
    """
    # Iterate once over the data and computes the error
    composition_loss = []
    for tensors in dataloader:
        _, _, losses = module(tensors, **kwargs)
        if isinstance(losses.composition_loss, dict):
            composition_reconstruction_loss = torch.stack(
                list(losses.composition_loss.values())
            ).sum(dim=0)
        else:
            composition_reconstruction_loss = losses.composition_loss
        composition_loss.append(composition_reconstruction_loss)

    composition_loss = torch.cat(composition_loss, dim=0)
    if return_mean:
        composition_loss = composition_loss.mean()
    return composition_loss


def compute_niche_error(
    module: Callable[[dict[str, Tensor | None], dict], tuple[Any, Any, NicheLossOutput]],
    dataloader: Iterator[dict[str, Tensor | None]],
    return_mean: bool = True,
    **kwargs,
) -> float | Tensor[float]:
    """Compute the niche state prediction error on the data.

    The  error is the negative log likelihood of the data (eta) given the latent
    variables.

    Parameters
    ----------
    module
        A callable (can be a :class:`~torch.nn.Module` instance) that takes a dictionary of
        :class:`~torch.Tensor`s as input and returns a tuple of three elements, where the last
        element is an instance of :class:`~scvi.external.scviva._module.NicheLossOutput`.
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
    The niche state prediction error of the data.
    """
    # Iterate once over the data and computes the error
    niche_loss = []
    for tensors in dataloader:
        _, _, losses = module(tensors, **kwargs)
        if isinstance(losses.niche_loss, dict):
            niche_reconstruction_loss = torch.stack(list(losses.niche_loss.values())).sum(dim=0)
        else:
            niche_reconstruction_loss = losses.niche_loss
        niche_loss.append(niche_reconstruction_loss)

    niche_loss = torch.cat(niche_loss, dim=0)
    if return_mean:
        niche_loss = niche_loss.mean()
    return niche_loss
