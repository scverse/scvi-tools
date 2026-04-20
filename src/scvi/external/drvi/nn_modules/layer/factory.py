from __future__ import annotations

from typing import TYPE_CHECKING

from scvi.external.drvi.nn_modules.layer.linear_layer import LinearLayer, StackedLinearLayer
from scvi.external.drvi.nn_modules.layer.structures import SimpleResidual

if TYPE_CHECKING:
    from typing import Any, Literal

    from torch import nn


class LayerFactory:
    """Abstract base class for creating neural network layers.

    This class provides a factory pattern for creating different types of
    neural network layers. It supports both normal layers and stacked layers
    with configurable architectures and residual connections.

    Parameters
    ----------
    intermediate_arch
        Architecture type for intermediate layers:
        - "SAME": Use the same architecture as defined in subclasses
        - "FC": Use fully connected layers
    residual_preferred
        Whether to wrap layers in residual connections when input and output
        dimensions match.

    Notes
    -----
    This is an abstract base class. Subclasses must implement:
    - `_get_normal_layer`: Creates normal layers
    - `_get_stacked_layer`: Creates stacked layers

    The factory pattern allows for easy switching between different layer
    architectures while maintaining a consistent interface.
    """

    def __init__(
        self, intermediate_arch: Literal["SAME", "FC"] = "SAME", residual_preferred: bool = False
    ) -> None:
        assert intermediate_arch in ["SAME", "FC"]

        self.intermediate_arch = intermediate_arch
        self.residual_preferred = residual_preferred

    def _get_normal_layer(
        self, d_in: int, d_out: int, bias: bool = True, **kwargs: Any
    ) -> nn.Module:
        """Create a normal layer (to be implemented by subclasses).

        Parameters
        ----------
        d_in
            Input dimension.
        d_out
            Output dimension.
        bias
            Whether to include bias term.
        **kwargs
            Additional keyword arguments.

        Returns
        -------
        nn.Module
            The created layer.

        Raises
        ------
        NotImplementedError
            This method must be implemented by subclasses.
        """
        raise NotImplementedError()

    def _get_stacked_layer(
        self, d_channel: int, d_in: int, d_out: int, bias: bool = True, **kwargs: Any
    ) -> nn.Module:
        """Create a stacked layer (to be implemented by subclasses).

        Parameters
        ----------
        d_channel
            Number of channels/splits.
        d_in
            Input dimension.
        d_out
            Output dimension.
        bias
            Whether to include bias term.
        **kwargs
            Additional keyword arguments.

        Returns
        -------
        nn.Module
            The created stacked layer.

        Raises
        ------
        NotImplementedError
            This method must be implemented by subclasses.
        """
        raise NotImplementedError()

    def get_normal_layer(
        self,
        d_in: int,
        d_out: int,
        bias: bool = True,
        intermediate_layer: bool | None = None,
        **kwargs: Any,
    ) -> nn.Module:
        """Create a normal layer with optional residual connection.

        Parameters
        ----------
        d_in
            Input dimension.
        d_out
            Output dimension.
        bias
            Whether to include bias term.
        intermediate_layer
            Whether this is an intermediate layer. If None, defaults to True.
        **kwargs
            Additional keyword arguments passed to layer creation.

        Returns
        -------
        nn.Module
            The created layer, optionally wrapped in a residual connection.

        Notes
        -----
        The layer creation logic depends on the `intermediate_arch` setting:
        - If `intermediate_layer=True` and `intermediate_arch="FC"`: Creates a linear layer
        - Otherwise: Uses the subclass-specific `_get_normal_layer` method

        If `residual_preferred=True` and `d_in == d_out`, the layer is wrapped
        in a `SimpleResidual` connection.

        Examples
        --------
        >>> factory = FCLayerFactory()
        >>> # Create intermediate layer with FC architecture
        >>> factory.intermediate_arch = "FC"
        >>> layer = factory.get_normal_layer(64, 128, intermediate_layer=True)
        >>> print(type(layer))  # <class 'torch.nn.modules.linear.Linear'>
        >>> # Create with residual connection
        >>> factory.residual_preferred = True
        >>> layer = factory.get_normal_layer(64, 64)
        >>> print(type(layer))  # <class 'drvi.nn_modules.layer.structures.SimpleResidual'>
        """
        if intermediate_layer is None:
            intermediate_layer = True
        if intermediate_layer and self.intermediate_arch == "FC":
            layer = LinearLayer(d_in, d_out, bias)
        elif (not intermediate_layer) or self.intermediate_arch == "SAME":
            layer = self._get_normal_layer(d_in, d_out, bias=True, **kwargs)
        else:
            raise NotImplementedError()

        if self.residual_preferred and (d_in == d_out):
            layer = SimpleResidual(layer)
        return layer

    def get_stacked_layer(
        self,
        d_channel: int,
        d_in: int,
        d_out: int,
        bias: bool = True,
        intermediate_layer: bool | None = None,
        **kwargs: Any,
    ) -> nn.Module:
        """Create a stacked layer with optional residual connection.

        Parameters
        ----------
        d_channel
            Number of channels/splits for the stacked layer.
        d_in
            Input dimension.
        d_out
            Output dimension.
        bias
            Whether to include bias term.
        intermediate_layer
            Whether this is an intermediate layer. If None, defaults to True.
        **kwargs
            Additional keyword arguments passed to layer creation.

        Returns
        -------
        nn.Module
            The created stacked layer, optionally wrapped in a residual connection.

        Notes
        -----
        The layer creation logic depends on the `intermediate_arch` setting:
        - If `intermediate_layer=True` and `intermediate_arch="FC"`: Creates a `StackedLinearLayer`
        - Otherwise: Uses the subclass-specific `_get_stacked_layer` method

        If `residual_preferred=True` and `d_in == d_out`, the layer is wrapped
        in a `SimpleResidual` connection.

        Stacked layers are useful for processing multiple splits of the input
        in parallel, which is common in disentanglement models.

        Examples
        --------
        >>> factory = FCLayerFactory()
        >>> # Create stacked layer with 4 splits
        >>> layer = factory.get_stacked_layer(4, 64, 128)
        >>> print(type(layer))  # <class 'drvi.nn_modules.layer.linear_layer.StackedLinearLayer'>
        >>> # Create with residual connection
        >>> factory.residual_preferred = True
        >>> layer = factory.get_stacked_layer(4, 64, 64)
        >>> print(type(layer))  # <class 'drvi.nn_modules.layer.structures.SimpleResidual'>
        """
        if intermediate_layer is None:
            intermediate_layer = True
        if intermediate_layer and self.intermediate_arch == "FC":
            layer = StackedLinearLayer(d_channel, d_in, d_out, bias)
        elif (not intermediate_layer) or self.intermediate_arch == "SAME":
            layer = self._get_stacked_layer(d_channel, d_in, d_out, bias=True, **kwargs)
        else:
            raise NotImplementedError()

        if self.residual_preferred and (d_in == d_out):
            layer = SimpleResidual(layer)
        return layer


class FCLayerFactory(LayerFactory):
    """Factory for creating fully connected neural network layers.

    This factory creates standard fully connected (linear) layers and stacked
    linear layers. It inherits from LayerFactory and implements the abstract
    methods to provide concrete layer creation functionality.

    Parameters
    ----------
    intermediate_arch
        Architecture type for intermediate layers:
        - "SAME": Use fully connected layers (same as "FC")
        - "FC": Use fully connected layers
    residual_preferred
        Whether to wrap layers in residual connections when input and output
        dimensions match.

    Notes
    -----
    This factory creates:
    - Normal layers: `LinearLayer` layers
    - Stacked layers: `StackedLinearLayer` for processing multiple splits

    The "SAME" and "FC" architectures are equivalent for this factory since
    both create fully connected layers.

    Examples
    --------
    >>> # Create factory with default settings
    >>> factory = FCLayerFactory()
    >>> # Create a simple linear layer
    >>> layer = factory.get_normal_layer(64, 128)
    >>> print(type(layer))  # <class 'torch.nn.modules.linear.Linear'>
    >>> # Create a stacked layer
    >>> stacked_layer = factory.get_stacked_layer(4, 64, 128)
    >>> print(
    ...     type(stacked_layer)
    ... )  # <class 'drvi.nn_modules.layer.linear_layer.StackedLinearLayer'>
    """

    def __init__(
        self, intermediate_arch: Literal["SAME", "FC"] = "SAME", residual_preferred: bool = False
    ) -> None:
        super().__init__(
            intermediate_arch=intermediate_arch, residual_preferred=residual_preferred
        )

    def _get_normal_layer(
        self, d_in: int, d_out: int, bias: bool = True, **kwargs: Any
    ) -> LinearLayer:
        """Create a fully connected layer.

        Parameters
        ----------
        d_in
            Input dimension.
        d_out
            Output dimension.
        bias
            Whether to include bias term.
        **kwargs
            Additional keyword arguments (ignored for linear layers).

        Returns
        -------
        LinearLayer
            A fully connected linear layer.

        Examples
        --------
        >>> factory = FCLayerFactory()
        >>> layer = factory._get_normal_layer(64, 128)
        >>> print(layer.weight.shape)  # torch.Size([128, 64])
        >>> print(layer.bias.shape)  # torch.Size([128])
        """
        return LinearLayer(d_in, d_out, bias=bias)

    def _get_stacked_layer(
        self, d_channel: int, d_in: int, d_out: int, bias: bool = True, **kwargs: Any
    ) -> StackedLinearLayer:
        """Create a stacked linear layer.

        Parameters
        ----------
        d_channel
            Number of channels/splits for the stacked layer.
        d_in
            Input dimension.
        d_out
            Output dimension.
        bias
            Whether to include bias term.
        **kwargs
            Additional keyword arguments (ignored for stacked linear layers).

        Returns
        -------
        StackedLinearLayer
            A stacked linear layer that processes multiple splits in parallel.

        Notes
        -----
        The `StackedLinearLayer` applies the same linear transformation to
        multiple splits of the input, which is useful for disentanglement
        models where different splits may represent different factors of variation.

        Examples
        --------
        >>> factory = FCLayerFactory()
        >>> layer = factory._get_stacked_layer(4, 64, 128)
        >>> print(layer.weight.shape)  # torch.Size([4, 128, 64])
        >>> print(layer.bias.shape)  # torch.Size([4, 128])
        """
        return StackedLinearLayer(d_channel, d_in, d_out, bias=bias)

    def __str__(self) -> str:
        """String representation of the factory.

        Returns
        -------
        str
            A string describing the factory configuration.
        """
        return f"FCLayerFactory(residual_preferred={self.residual_preferred})"
