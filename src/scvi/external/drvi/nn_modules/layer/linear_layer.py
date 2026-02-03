from __future__ import annotations

import math
from typing import TYPE_CHECKING

import torch
from torch import nn
from torch.nn import functional as F

if TYPE_CHECKING:
    from typing import Any


class LinearLayer(nn.Linear):
    def forward(self, x: torch.Tensor, output_subset: torch.Tensor | None = None) -> torch.Tensor:
        if output_subset is None:
            # x: (..., i) -> output: (..., o)
            return super().forward(x)
        elif output_subset.dim() == 1:
            # x: (..., i) -> output_subset: (o_subset)
            bias = self.bias[output_subset] if self.bias is not None else None  # (o_subset)
            weight = self.weight[output_subset]  # (o_subset, i)
            return F.linear(x, weight, bias)  # (..., i) -> (..., o_subset)
        else:
            raise NotImplementedError()


class StackedLinearLayer(nn.Module):
    """A stacked linear layer that applies multiple linear transformations in parallel.

    This layer applies a linear transformation to multiple channels/splits
    of the input. It's particularly useful in additive decoders where
    different splits should be calculated in parallel.

    Parameters
    ----------
    n_channels
        Number of channels/splits to process in parallel.
    in_features
        Number of input features per channel.
    out_features
        Number of output features per channel.
    bias
        Whether to include bias terms for each channel.
    device
        Device to place the layer on.
    dtype
        Data type for the layer parameters.

    Notes
    -----
    The layer maintains separate weight and bias parameters for each channel:
    - Weight shape: (n_channels, in_features, out_features)
    - Bias shape: (n_channels, out_features) if bias=True, None otherwise

    The forward pass applies the transformation to each channel independently:
    output[b, c, o] = sum_i(x[b, c, i] * weight[c, i, o]) + bias[c, o]

    This is equivalent to applying n_channels separate linear layers in parallel,
    which is more efficient than using separate nn.Linear layers.

    Examples
    --------
    >>> import torch
    >>> # Create a stacked linear layer with 4 channels
    >>> layer = StackedLinearLayer(n_channels=4, in_features=64, out_features=128)
    >>> # Input shape: (batch_size, n_channels, in_features)
    >>> x = torch.randn(32, 4, 64)
    >>> # Forward pass
    >>> output = layer(x)
    >>> print(output.shape)  # torch.Size([32, 4, 128])
    >>> # Each channel has its own parameters
    >>> print(layer.weight.shape)  # torch.Size([4, 64, 128])
    >>> print(layer.bias.shape)  # torch.Size([4, 128])
    """

    __constants__ = ["n_channels", "in_features", "out_features"]
    n_channels: int
    in_features: int
    out_features: int
    weight: torch.Tensor
    bias: torch.Tensor | None

    def __init__(
        self,
        n_channels: int,
        in_features: int,
        out_features: int,
        bias: bool = True,
        device: Any = None,
        dtype: Any = None,
    ) -> None:
        factory_kwargs = {"device": device, "dtype": dtype}
        super().__init__()
        self.n_channels = n_channels
        self.in_features = in_features
        self.out_features = out_features
        self.weight = nn.Parameter(torch.empty((n_channels, in_features, out_features), **factory_kwargs))
        if bias:
            self.bias = nn.Parameter(torch.empty(n_channels, out_features, **factory_kwargs))
        else:
            self.register_parameter("bias", None)
        self.reset_parameters()

    def reset_parameters(self) -> None:
        """Reset the layer parameters to their initial values.

        This method reinitializes both weights and biases using the same
        initialization strategy as the default nn.Linear layer.

        Notes
        -----
        The initialization follows PyTorch's default linear layer initialization:
        - Weights: Uniform distribution in [-1/sqrt(in_features), 1/sqrt(in_features)]
        - Biases: Uniform distribution in [-1/sqrt(in_features), 1/sqrt(in_features)]

        This ensures that the variance of the output is approximately preserved
        across the layer.
        """
        self._init_weight()
        self._init_bias()

    def _init_weight(self) -> None:
        """Initialize the weight parameters.

        Notes
        -----
        Uses the same initialization as default nn.Linear:
        Uniform distribution in [-1/sqrt(in_features), 1/sqrt(in_features)]

        This initialization helps maintain the variance of activations
        across the network, which is important for training stability.
        """
        # Same as default nn.Linear (https://github.com/pytorch/pytorch/issues/57109)
        fan_in = self.in_features
        bound = 1 / math.sqrt(fan_in)
        nn.init.uniform_(self.weight, -bound, bound)

    def _init_bias(self) -> None:
        """Initialize the bias parameters.

        Notes
        -----
        Uses the same initialization as default nn.Linear:
        Uniform distribution in [-1/sqrt(in_features), 1/sqrt(in_features)]

        The bias initialization is independent of the weight initialization
        and helps ensure that the layer can learn appropriate offsets.
        """
        if self.bias is not None:
            fan_in = self.in_features
            bound = 1 / math.sqrt(fan_in)
            nn.init.uniform_(self.bias, -bound, bound)

    def forward(self, x: torch.Tensor, output_subset: torch.Tensor | None = None) -> torch.Tensor:
        r"""Forward pass through the stacked linear layer.

        Parameters
        ----------
        x
            Input tensor with shape (batch_size, n_channels, in_features).
        output_subset
            Subset of outputs to provide in the output.

        Returns
        -------
        torch.Tensor
            Output tensor with shape (batch_size, n_channels, out_features).

        Notes
        -----
        The forward pass applies the linear transformation to each channel:

        .. math::
            \text{output}[b, c, o] = \\sum_{i} \text{input}[b, c, i] \\cdot \text{weight}[c, i, o] + \text{bias}[c, o]

        where:
        - b: batch index
        - c: channel index
        - i: input feature index
        - o: output feature index

        The computation is performed efficiently using torch.einsum, which
        is equivalent to applying n_channels separate linear transformations
        in parallel.

        Examples
        --------
        >>> import torch
        >>> # Create layer
        >>> layer = StackedLinearLayer(n_channels=3, in_features=10, out_features=5)
        >>> # Input: batch_size=2, n_channels=3, in_features=10
        >>> x = torch.randn(2, 3, 10)
        >>> # Forward pass
        >>> output = layer(x)
        >>> print(output.shape)  # torch.Size([2, 3, 5])
        """
        if True:
            if output_subset is None or output_subset.dim() == 1:
                # weight: (c, i, o), bias: (c, o)
                # x: (b, c, i), output_subset: (o_subset) -> output: (b, c, o_subset)
                weight = self.weight if output_subset is None else self.weight[:, :, output_subset]  # (c, i, o_subset)
                # slower: mm = torch.einsum("bci,cio->bco", x, weight)
                mm = torch.bmm(x.transpose(0, 1), weight).transpose(0, 1)  # (b, c, o_subset)
                if self.bias is not None:
                    bias = self.bias if output_subset is None else self.bias[:, output_subset]  # (c, o_subset)
                    mm = mm + bias  # They (bco, co) will broadcast well
                return mm
            else:
                raise NotImplementedError()

    def extra_repr(self) -> str:
        """String representation for printing the layer.

        Returns
        -------
        str
            A string containing the layer's configuration parameters.

        Notes
        -----
        This method is used by PyTorch's __repr__ method to provide
        a detailed string representation of the layer, which is useful
        for debugging and understanding the layer's configuration.

        Examples
        --------
        >>> layer = StackedLinearLayer(n_channels=4, in_features=64, out_features=128, bias=True)
        >>> print(layer)
        StackedLinearLayer(in_features=64, out_features=128, n_channels=4, bias=True)
        """
        return (
            f"in_features={self.in_features}, out_features={self.out_features}, "
            f"n_channels={self.n_channels}, bias={self.bias is not None}"
        )
