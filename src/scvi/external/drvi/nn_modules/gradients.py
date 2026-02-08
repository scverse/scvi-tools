import torch


class GradScale(torch.autograd.Function):
    """Custom autograd function for scaling gradients during backward pass."""

    @staticmethod
    def forward(ctx, x, scale):
        """Forward pass that stores the scale factor.

        Parameters
        ----------
        ctx
            Context object for storing information between forward and backward.
        x
            Input tensor (unchanged in forward pass).
        scale
            Scale factor to apply during backward pass.

        Returns
        -------
        torch.Tensor
            Input tensor unchanged.
        """
        ctx.scale = scale
        return x  # forward pass unchanged

    @staticmethod
    def backward(ctx, grad_output):
        """Backward pass that scales the gradient.

        Parameters
        ----------
        ctx
            Context object containing the scale factor.
        grad_output
            Gradient from the next layer.

        Returns
        -------
        tuple
            Scaled gradient and None (no gradient for scale parameter).
        """
        return grad_output * ctx.scale, None  # scale gradient only


def grad_scale(x, scale):
    """Apply gradient scaling to a tensor.

    Parameters
    ----------
    x
        Input tensor.
    scale
        Scale factor to apply to gradients during backward pass.

    Returns
    -------
    torch.Tensor
        Tensor with gradient scaling applied.
    """
    return GradScale.apply(x, scale)


class GradientScaler(torch.nn.Module):
    """Module that scales gradients during backward pass.

    This module wraps the grad_scale function to provide a PyTorch module
    interface for gradient scaling. Useful for controlling gradient magnitudes
    in specific parts of the network.

    Parameters
    ----------
    scale
        Scale factor to apply to gradients during backward pass.
    """

    def __init__(self, scale: float):
        """Initialize the gradient scaler.

        Parameters
        ----------
        scale
            Scale factor to apply to gradients during backward pass.
        """
        super().__init__()
        self.register_buffer("scale", torch.tensor(scale, dtype=torch.float32))

    def forward(self, x):
        """Forward pass that applies gradient scaling.

        Parameters
        ----------
        x
            Input tensor.

        Returns
        -------
        torch.Tensor
            Input tensor with gradient scaling registered.
        """
        return grad_scale(x, self.scale)

    def extra_repr(self):
        """String representation of the module.

        Returns
        -------
        str
            String showing the scale factor.
        """
        return f"scale={self.scale.item()}"
