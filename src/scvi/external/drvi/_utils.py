from __future__ import annotations

import math
from typing import TYPE_CHECKING

import torch
from torch import nn

if TYPE_CHECKING:
    from typing import Any


class StackedLinearLayer(nn.Module):
    """A parallel stacked linear layer that applies multiple linear transformations in parallel.

    This layer applies a linear transformation to multiple stacks/splits
    of the input. It's particularly useful in additive decoders where
    different splits should be calculated in parallel.

    Parameters
    ----------
    n_stacks
        Number of stacks/splits to process in parallel.
    in_features
        Number of input features per stack.
    out_features
        Number of output features per stack.
    bias
        Whether to include bias terms for each stack.
    device
        Device to place the layer on.
    dtype
        Data type for the layer parameters.

    Notes
    -----
    The layer maintains separate weight and bias parameters for each stack:
    - Weight shape: (n_stacks, out_features, in_features)
    - Bias shape: (n_stacks, out_features) if bias=True, None otherwise

    The forward pass applies the transformation to each stack independently:
    output[b, s, o] = sum_i(x[b, s, i] * weight[s, o, i]) + bias[s, o]

    This is equivalent to applying n_stacks separate linear layers in parallel,
    which is more efficient than using separate nn.Linear layers.

    Examples
    --------
    >>> import torch
    >>> # Create a stacked linear layer with 4 stacks
    >>> layer = StackedLinearLayer(n_stacks=4, in_features=64, out_features=128)
    >>> # Input shape: (batch_size, n_stacks, in_features)
    >>> x = torch.randn(32, 4, 64)
    >>> # Forward pass
    >>> output = layer(x)
    >>> print(output.shape)  # torch.Size([32, 4, 128])
    >>> # Each stack has its own parameters
    >>> print(layer.weight.shape)  # torch.Size([4, 128, 64])
    >>> print(layer.bias.shape)  # torch.Size([4, 128])
    """

    __constants__ = ["n_stacks", "in_features", "out_features"]
    n_stacks: int
    in_features: int
    out_features: int
    weight: torch.Tensor
    bias: torch.Tensor | None

    def __init__(
        self,
        n_stacks: int,
        in_features: int,
        out_features: int,
        bias: bool = True,
        device: Any = None,
        dtype: Any = None,
    ) -> None:
        factory_kwargs = {"device": device, "dtype": dtype}
        super().__init__()
        self.n_stacks = n_stacks
        self.in_features = in_features
        self.out_features = out_features
        self.weight = nn.Parameter(
            torch.empty((n_stacks, out_features, in_features), **factory_kwargs)
        )
        if bias:
            self.bias = nn.Parameter(torch.empty(n_stacks, out_features, **factory_kwargs))
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

    def forward(
        self,
        x: torch.Tensor,
        output_subset: torch.Tensor | None = None,
        stack_subset: torch.Tensor | None = None,
    ) -> torch.Tensor:
        r"""Forward pass through the stacked linear layer.

        Parameters
        ----------
        x
            Input tensor with shape (..., n_stacks, in_features). Any number of
            leading (batch) dimensions is supported; the transformation is only
            applied to the last two dimensions (stacks and features).
        output_subset
            Subset of outputs to provide in the output.
        stack_subset
            Indices for stacks in operation.

        Returns
        -------
        torch.Tensor
            Output tensor with shape (..., n_stacks, out_features), matching the
            leading dimensions of the input.

        Notes
        -----
        The forward pass applies the linear transformation to each stack:

        where:
        - ...: arbitrary leading (batch) dimensions
        - s: stack index
        - i: input feature index
        - o: output feature index

        The computation is performed efficiently using torch.bmm or broadcasting.

        Examples
        --------
        >>> import torch
        >>> # Create layer
        >>> layer = StackedLinearLayer(n_stacks=3, in_features=10, out_features=5)
        >>> # Input: batch_size=2, n_stacks=3, in_features=10
        >>> x = torch.randn(2, 3, 10)
        >>> # Forward pass
        >>> output = layer(x)
        >>> print(output.shape)  # torch.Size([2, 3, 5])
        >>> # Extra leading dimensions are also supported
        >>> x = torch.randn(7, 2, 3, 10)
        >>> layer(x).shape  # torch.Size([7, 2, 3, 5])
        """
        if stack_subset is None:
            if output_subset is None or output_subset.dim() == 1:
                # weight: (s, o, i), bias: (s, o)
                # x: (..., s, i), output_subset: (o_subset) -> output: (..., s, o_subset)
                weight = (
                    self.weight if output_subset is None else self.weight[:, output_subset]
                )  # (s, o_subset, i)
                # Flatten the leading dims into a single batch so we can use the fast
                # batched matmul (bmm), then restore the original leading dims.
                *lead_shape, n_stacks, in_features = x.shape
                x_flat = x.reshape(-1, n_stacks, in_features)  # (b, s, i)
                # slower: mm = torch.einsum("bsi,soi->bso", x_flat, weight)
                mm = torch.bmm(x_flat.transpose(0, 1), weight.transpose(1, 2)).transpose(
                    0, 1
                )  # (b, s, o_subset)
                mm = mm.reshape(*lead_shape, n_stacks, mm.shape[-1])  # (..., s, o_subset)
                if self.bias is not None:
                    bias = (
                        self.bias if output_subset is None else self.bias[:, output_subset]
                    )  # (s, o_subset)
                    mm = mm + bias  # (..., s, o_subset) and (s, o_subset) will broadcast well
                return mm
            else:
                raise NotImplementedError()
        else:
            # stack_subset: (..., s_subset)
            # x: (..., s_subset, i), output_subset: (o_subset) -> output: (..., s_subset, o_subset)
            weight = self.weight[stack_subset]  # (..., s_subset, o, i)
            bias = self.bias[stack_subset] if self.bias is not None else None  # (..., s_subset, o)

            if output_subset is None:
                pass
            elif output_subset.dim() == 1:
                weight = weight[..., output_subset, :]  # (..., s_subset, o_subset, i)
                bias = (
                    bias[..., output_subset] if bias is not None else None
                )  # (..., s_subset, o_subset)
            else:
                raise NotImplementedError
            mm = torch.matmul(x.unsqueeze(-2), weight.transpose(-1, -2)).squeeze(
                -2
            )  # (..., s_subset, o_subset)
            if bias is not None:
                mm = mm + bias  # (..., s_subset, o_subset)
            return mm

    def extra_repr(self) -> str:
        """String representation for printing the layer."""
        return (
            f"in_features={self.in_features}, out_features={self.out_features}, "
            f"n_stacks={self.n_stacks}, bias={self.bias is not None}"
        )
