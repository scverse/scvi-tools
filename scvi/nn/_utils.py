from typing import Optional

import jax.numpy as jnp
import torch

from scvi._types import Tensor

from .jax._utils import _one_hot_jax
from .torch._utils import _one_hot_torch


def one_hot(indexes: Tensor, n_classes: Optional[int] = None) -> Tensor:
    """
    One-hot encode a tensor of category indexes.

    Accepts either a :class:`~torch.Tensor` or a :class:`~jax.numpy.ndarray`.

    Parameters
    ----------
    indexes
        A :class:`~scvi._types.Tensor` of shape `(n_samples,)` or `(n_samples, 1)`
        containing the category index for each sample.
    n_classes
        The number of categories. If `None`, the number of categories is inferred from
        the maximum value in `indexes`.

    Returns
    -------
    oh
        A :class:`~scvi._types.Tensor` of shape `(n_samples, n_classes)` containing the
        one-hot encoding of `indexes`.
    """
    if isinstance(indexes, torch.Tensor):
        oh = _one_hot_torch(indexes, n_classes=n_classes)
    elif isinstance(indexes, jnp.ndarray):
        oh = _one_hot_jax(indexes, n_classes=n_classes)
    else:
        raise TypeError("`indexes` must be a `torch.Tensor` or `jax.numpy.ndarray`")
    return oh
