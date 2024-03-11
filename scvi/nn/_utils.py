import warnings

import torch

from scvi import settings


def _one_hot(index: torch.Tensor, n_cat: int) -> torch.Tensor:
    """Wrapper for PyTorch's one-hot encoding.

    Assumes an input of dimensions N x ... x 1, and collapses last
    dimension to give output dimensions N x ... x n_cat.
    """
    if index.shape[-1] != 1:
        raise ValueError(f"Expected input with dimensions [N, ..., 1]. Received {index.shape}")

    onehot = torch.nn.functional.one_hot(index, n_cat).squeeze(-2)
    return onehot.type(torch.float32)


def one_hot(index: torch.Tensor, n_cat: int) -> torch.Tensor:
    """One hot a tensor of categories."""
    # TODO: remove in 1.3.0
    warnings.warn(
        "The `one_hot` function is deprecated and will be removed in scvi-tools 1.3. "
        "Please use the `one_hot` function in PyTorch instead.",
        DeprecationWarning,
        stacklevel=settings.warnings_stacklevel,
    )
    onehot = torch.zeros(index.size(0), n_cat, device=index.device)
    onehot.scatter_(1, index.type(torch.long), 1)
    return onehot.type(torch.float32)
