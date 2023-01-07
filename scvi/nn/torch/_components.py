from typing import Optional

import torch
from chex import dataclass
from torch import nn


@dataclass
class Module(nn.Module):
    """Dummy class to allow PyTorch modules to be dataclasses."""

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.setup()

    def setup(self):
        """Setup the module."""


class BatchNorm(nn.Module):
    """Batch normalization layer."""

    def __init__(self, *args, **kwargs):
        super().__init__()
        self._module = nn.BatchNorm1d(*args, **kwargs)

    def forward(self, x: torch.Tensor) -> torch.Tensor:
        """
        Forward pass through the module.

        Parameters
        ----------
        x
            Input :class:`~torch.Tensor` of shape `(n_samples, n_input)` or
            `(n_layers, n_samples, n_input)`. If 3-dimensional, the operation
            is applied to each layer separately and then concatenated.

        Returns
        -------
        out
            Output :class:`~torch.Tensor` of shape `(n_samples, n_input)` or
            `(n_layers, n_samples, n_input)`.
        """

        def _forward(x_: torch.Tensor) -> torch.Tensor:
            # (n, d) -> (1, n, d)
            return self._module(x_).unsqueeze(0)

        if x.dim() == 3:
            out = torch.cat(map(_forward, x), dim=0)
        else:
            out = self._module(x)
        return out

    def __getattr__(self, name):
        return getattr(self._module, name)


class Linear(nn.Module):
    """
    Linear layer that allows for embeddings to be injected.

    Parameters
    ----------
    args
        Arguments for :class:`~torch.nn.Linear`.
    inject_embedding
        Whether to inject embeddings.
    **kwargs
        Keyword arguments for :class:`~torch.nn.Linear`.
    """

    def __init__(self, *args, inject_embedding: bool = False, **kwargs):
        super().__init__()
        self._inject_embedding = inject_embedding
        self._module = nn.Linear(*args, **kwargs)

    def forward(
        self, x: torch.Tensor, embedding: Optional[torch.Tensor] = None
    ) -> torch.Tensor:
        """
        Forward pass through the module.

        Parameters
        ----------
        x
            Input :class:`~torch.Tensor` of shape `(n_samples, n_input)` or
            `(n_layers, n_samples, n_input)`.
        embedding
            Embedding :class:`~torch.Tensor` of shape `(n_samples, embedding_dim)`.

        Returns
        -------
        out
            Output :class:`~torch.Tensor` of shape `(n_samples, n_output)` or
            `(n_layers, n_samples, n_output)`.
        """
        if not self._inject_embedding and embedding is not None:
            raise ValueError("Embedding not allowed in this layer.")

        if self._inject_embedding and embedding is not None:
            if x.dim() == 3:
                # (n_samples, embedding_dim) -> (n_layers, n_samples, embedding_dim)
                embedding = embedding.unsqueeze(0).expand(x.size(0), -1, -1)
            x = torch.cat((x, embedding), dim=-1)

        out = self._module(x)
        return out

    def __getattr__(self, name):
        return getattr(self._module, name)
