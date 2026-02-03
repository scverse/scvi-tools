from __future__ import annotations

import torch
from torch import nn


class SimpleResidual(nn.Module):
    """Simple residual connection that adds the input to the layer output.

    This module implements a basic residual connection by adding the input
    directly to the output of the wrapped layer. It's useful for creating
    skip connections that help with gradient flow in deep networks.

    Parameters
    ----------
    layer : nn.Module
        The layer to be wrapped with a residual connection.

    Notes
    -----
    The forward pass computes: output = input + layer(input)

    This is the simplest form of residual connection, where the input is
    added directly to the layer's output. For this to work properly, the
    layer should preserve the input shape (i.e., input.shape == layer(input).shape).

    Examples
    --------
    >>> import torch
    >>> from torch import nn
    >>> # Create a simple linear layer
    >>> linear_layer = nn.Linear(64, 64)
    >>> # Wrap it with a residual connection
    >>> residual_layer = SimpleResidual(linear_layer)
    """

    def __init__(self, layer: nn.Module):
        """Initialize the residual connection.

        Parameters
        ----------
        layer
            The layer to be wrapped with a residual connection.

        Notes
        -----
        The layer should preserve the input shape for the residual connection
        to work properly. This means layer(input).shape should equal input.shape.
        """
        super().__init__()
        self.layer = layer

    def forward(self, x: torch.Tensor) -> torch.Tensor:
        r"""Forward pass through the residual connection.

        Parameters
        ----------
        x
            Input tensor. The shape should be compatible with the wrapped layer.

        Returns
        -------
        torch.Tensor
            Output tensor with the same shape as the input.

        Notes
        -----
        The forward pass computes the residual connection:

        .. math::
            \text{output} = \text{input} + \text{layer}(\text{input})

        This is equivalent to adding a skip connection from the input to the
        output of the wrapped layer.
        """
        return x + self.layer(x)
