from typing import List, Literal, Optional

import flax
import torch

from scvi.utils import InvalidParameterError

from .jax._fc_layers import JaxFCLayers
from .torch._fc_layers import TorchFCLayers


class FCLayers:
    """
    A helper class to build fully-connected layers for scvi-tools models.

    Use the :meth:`~scvi.nn.FCLayers.setup` method to instantiate the correct backend.
    """

    @staticmethod
    def setup(
        *,
        backend: Literal["torch", "jax"] = "torch",
        n_input: Optional[int] = None,
        n_output: Optional[int] = None,
        n_layers: Optional[int] = None,
        n_hidden: Optional[int] = None,
        layers_dim: Optional[List[int]] = None,
        bias: bool = True,
        dropout_rate: Optional[float] = 0.1,
        norm: Optional[Literal["batch", "layer", "group"]] = "batch",
        norm_kwargs: Optional[dict] = None,
        activation: Optional[str] = "relu",
        activation_kwargs: Optional[dict] = None,
        residual: bool = False,
        n_cat_list: Optional[List[int]] = None,
        inject_covariates: bool = True,
        training: Optional[bool] = None,
    ):
        """
        Helper method to instantiate a fully-connected network.

        Validates inputs before passing them to the correct backend.

        Parameters
        ----------
        backend
            The backend to use. One of the following:

            * ``'torch'``: :class:`~scvi.nn.TorchFCLayers`
            * ``'jax'``: :class:`~scvi.nn.JaxFCLayers`
        n_input
            The number of input features. Ignored if `backend` is `'jax'`.
        n_output
            The number of output features. Must be provided if `layers_dim` is not.
        n_layers
            The number of hidden layers including the output layer. Must be provided if
            `layers_dim` is not.
        n_hidden
            The number of nodes per hidden layer. Must be provided if `layers_dim` is not.
        layers_dim
            If provided, overrides `n_output`, `n_layers`, and `n_hidden`. Values correspond
            to the dimensions of each fully-connected layer such that the last value is the
            output dimension.
        bias
            Whether to include bias terms in the layers.
        dropout_rate
            The dropout rate to apply to the layers. If `None`, no dropout is applied.
        norm
            The normalization to apply to the layers. One of the following:

            * ``'batch'``: :class:`~torch.nn.BatchNorm1d` or :class:`~jax.nn.BatchNorm`
            * ``'layer'``: :class:`~torch.nn.LayerNorm` or :class:`~jax.nn.LayerNorm`
            * ``'group'``: :class:`~torch.nn.GroupNorm` or :class:`~jax.nn.GroupNorm`
            * ``None``: no normalization
        norm_kwargs
            Keyword arguments to pass to the normalization layer.
        activation
            The activation to apply to the layers. Must correspond to a function in
            :attr:`torch.nn.functional` or :attr:`jax.nn`. If `None`, no activation is used.
        activation_kwargs
            Keyword arguments to pass to the activation layer.
        residual
            Whether to include residual connections between hidden layers. Cannot be used if
            `layers_dim` is provided since uniform hidden dimensions cannot be guaranteed.
        n_cat_list
            A list of integers specifying the number of categories for each covariate. Each
            category will be included using a one-hot encoding.
        inject_covariates
            Whether to inject covariates into each hidden layer in addition to the input
            layer.
        training
            Whether the module is in training mode.

        Returns
        -------
        :class:`~scvi.nn.TorchFCLayers` or :class:`~scvi.nn.JaxFCLayers`
        """
        valid_backends = ["torch", "jax"]
        if backend not in valid_backends:
            raise InvalidParameterError(
                param="backend", value=backend, valid=valid_backends
            )

        # validate layer arguments
        if not all([n_input, n_output, n_layers, n_hidden]) and not layers_dim:
            raise ValueError(
                "Must specify either (`n_input`, `n_output`, `n_layers`, `n_hidden`) "
                "or `layers_dim`."
            )
        if layers_dim and residual:
            raise ValueError("Cannot set `residual=True` if `layers_dim` is provided.")
        if layers_dim is None:
            layers_dim = [n_hidden] * (n_layers - 1) + [n_output]
            if backend == "torch":
                # jax backend infers input dimension from data
                layers_dim = [n_input] + layers_dim

        # validate normalization arguments
        valid_norms = ["batch", "layer", "group", None]
        if norm not in valid_norms:
            raise InvalidParameterError(param="norm", value=norm, valid=valid_norms)
        norm_kwargs = norm_kwargs or {}

        # validate activation arguments
        if (
            activation
            and backend == "torch"
            and activation not in dir(torch.nn.functional)
        ):
            raise InvalidParameterError(
                param="activation",
                value=activation,
                additional_message="Must be a function in `torch.nn.functional`.",
            )
        elif activation and backend == "jax" and activation not in dir(flax.linen):
            raise InvalidParameterError(
                param="activation",
                value=activation,
                additional_message="Must be a function in `flax.linen`.",
            )
        activation_kwargs = activation_kwargs or {}

        # ignore covariates with only one category
        n_cat_list = n_cat_list or []
        n_cat_list = [n_cat if n_cat > 1 else 0 for n_cat in n_cat_list]

        kwargs = {
            "layers_dim": layers_dim,
            "bias": bias,
            "dropout_rate": dropout_rate,
            "norm": norm,
            "norm_kwargs": norm_kwargs,
            "activation": activation,
            "activation_kwargs": activation_kwargs,
            "residual": residual,
            "n_cat_list": n_cat_list,
            "inject_covariates": inject_covariates,
            "training": training,
        }
        return TorchFCLayers(**kwargs) if backend == "torch" else JaxFCLayers(**kwargs)
