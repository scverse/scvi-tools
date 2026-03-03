from __future__ import annotations

from typing import TYPE_CHECKING

from torch import nn

if TYPE_CHECKING:
    from typing import Any


def freezable(base_norm_class):
    """Decorator to create freezable versions of normalization layers.

    This decorator creates a wrapper class around a base normalization class
    that adds freezing capability. When frozen, the normalization layer
    operates in evaluation mode regardless of the model's training state.

    Parameters
    ----------
    base_norm_class
        The base normalization class to make freezable.

    Returns
    -------
    type
        A new class that inherits from base_norm_class with freezing capability.

    Notes
    -----
    The freezable wrapper adds:
    - A `_freeze` attribute to track freeze status
    - A `freeze()` method to control freezing
    - Modified `forward()` method that respects freeze status

    When frozen, the normalization layer:
    - Temporarily switches to evaluation mode
    - Uses running statistics instead of batch statistics
    - Prevents parameter updates

    Examples
    --------
    >>> import torch
    >>> from torch import nn
    >>> # Create a freezable batch norm
    >>> FreezableBN = freezable(nn.BatchNorm1d)
    >>> bn = FreezableBN(10)
    >>> # Test normal operation
    >>> x = torch.randn(5, 10)
    >>> output1 = bn(x)
    >>> # Freeze the layer
    >>> bn.freeze(True)
    >>> output2 = bn(x)
    >>> # The outputs will be different due to different normalization behavior
    """

    class FreezableNormClass(base_norm_class):
        """Freezable wrapper for normalization layers.

        This class adds freezing capability to normalization layers by
        temporarily switching to evaluation mode when frozen.

        Parameters
        ----------
        *args
            Arguments passed to the base normalization class.
        **kwargs
            Keyword arguments passed to the base normalization class.
        """

        def __init__(self, *args: Any, **kwargs: Any) -> None:
            self._freeze = False
            super().__init__(*args, **kwargs)

        def freeze(self, freeze_status: bool = True) -> None:
            """Set the freeze status of the normalization layer.

            Parameters
            ----------
            freeze_status
                Whether to freeze the normalization layer.

            Notes
            -----
            When frozen, the layer will operate in evaluation mode during
            forward passes, using running statistics instead of batch statistics.
            This is useful for inference or when you want to prevent the
            normalization parameters from being updated.
            """
            self._freeze = freeze_status

        def forward(self, *args: Any, **kwargs: Any) -> Any:
            """Forward pass with freeze-aware behavior.

            Parameters
            ----------
            *args
                Arguments passed to the base normalization forward method.
            **kwargs
                Keyword arguments passed to the base normalization forward method.

            Returns
            -------
            torch.Tensor
                Output of the normalization layer.

            Notes
            -----
            If the layer is frozen, it temporarily switches to evaluation mode
            for the forward pass, then restores the original training state.
            This ensures that frozen layers use running statistics regardless
            of the model's overall training state.
            """
            training_status = self.training
            if self._freeze:
                self.train(False)
            result = super().forward(*args, **kwargs)
            self.train(training_status)
            return result

    return FreezableNormClass


FreezableBatchNorm1d = freezable(nn.BatchNorm1d)
FreezableLayerNorm = freezable(nn.LayerNorm)
