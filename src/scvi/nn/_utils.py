import warnings

import torch
from torch import nn

from scvi import settings


def one_hot(index: torch.Tensor, n_cat: int) -> torch.Tensor:
    """One hot a tensor of categories."""
    # TODO: remove in 1.3.0
    warnings.warn(
        "The `one_hot` function is deprecated in v1.2 and will be removed in v1.3. "
        "Please use the `one_hot` function in PyTorch instead.",
        DeprecationWarning,
        stacklevel=settings.warnings_stacklevel,
    )
    onehot = torch.zeros(index.size(0), n_cat, device=index.device)
    onehot.scatter_(1, index.type(torch.long), 1)
    return onehot.type(torch.float32)


class ExpActivation(nn.Module):
    def __init__(self):
        super().__init__()

    def forward(self, input):
        return torch.exp(input)
