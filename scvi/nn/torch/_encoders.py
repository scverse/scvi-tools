from typing import List, Literal, Optional, Tuple

import torch
from torch import nn
from torch.nn import functional as F

from ._components import Module
from ._fc_layers import TorchFCLayers


class TorchEncoder(Module):
    """
    Encodes data of `n_input` dimensions into `n_latent` dimensions.

    Parameters
    ----------
    n_input
        The number of input features. Must be provided.
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
        The activation to apply to the layers. If `None`, no activation is used.
    activation_kwargs
        Keyword arguments to pass to the activation layer.
    residual
        Whether to include residual connections between hidden layers. Cannot be used if
        `layers_dim` is provided since uniform hidden dimensions cannot be guaranteed.
    n_classes_list
        A list of integers specifying the number of classes for each covariate. Each
        class will be included using a one-hot encoding.
    embedding_dim
        The dimension of embeddings to use for encoding covariates. If `None`, a
        one-hot encoding is used.
    inject_covariates
        Whether to inject covariates into each hidden layer in addition to the input
        layer.
    training
        Whether the module is in training mode.
    distribution
        The distribution to use for the output. One of the following:

        * ``'normal'``: Normal distribution
        * ``'ln'``: Log-normal distribution
    variance_eps
        Minumum variance for numerical stability.
    variance_activation
        The activation to apply to the variance to ensure positivity. Must correspond to
        a function in :attr:`torch.nn.functional`.
    variance_activation_kwargs
        Keyword arguments to pass to the variance activation.
    return_distribution
        Whether to return the :class:`~torch.distributions.Distribution` object instead
        of its parameters.
    """

    n_input: int
    layers_dim: List[int]
    n_latent: int
    bias: bool
    dropout_rate: float
    norm: Literal["batch", "layer", "group", None]
    norm_kwargs: dict
    activation: str
    activation_kwargs: dict
    residual: bool
    n_classes_list: List[int]
    embedding_dim: Optional[int]
    inject_covariates: bool
    training: Optional[bool]
    distribution: Literal["normal", "ln"]
    variance_eps: float
    variance_activation: str
    variance_activation_kwargs: dict
    return_distribution: bool

    def setup(self):
        """Setup."""
        _kwargs = {
            "n_input": self.n_input,
            "layers_dim": self.layers_dim,
            "bias": self.bias,
            "dropout_rate": self.dropout_rate,
            "norm": self.norm,
            "norm_kwargs": self.norm_kwargs,
            "activation": self.activation,
            "activation_kwargs": self.activation_kwargs,
            "residual": self.residual,
            "n_classes_list": self.n_classes_list,
            "embedding_dim": self.embedding_dim,
            "inject_covariates": self.inject_covariates,
            "training": self.training,
        }
        self._module = TorchFCLayers(**_kwargs)
        self._mean = nn.Linear(self.layers_dim[-1], self.n_latent)
        self._variance = nn.Linear(self.layers_dim[-1], self.n_latent)

        self.z_transformation = lambda x: x
        if self.distribution == "ln":
            self.z_transformation = F.softmax

    def forward(
        self, x: torch.Tensor, indexes_list: Optional[List[torch.Tensor]] = None
    ) -> Tuple[torch.Tensor, Optional[torch.distribution.Distribution]]:
        """
        Forward pass through the module.

        Parameters
        ----------
        x
            Input :class:`~torch.Tensor` of shape `(n_samples, n_input)` or
            `(n_layers, n_samples, n_input)`.
        indexes_list
            List of :class:`~torch.Tensor` of shape `(n_samples, 1)` or
            `(n_samples, n_classes)` containing the categorical membership of each
            sample for a given covariate.

        Returns
        -------
        Tuple
        """
        q = self._module(x, indexes_list=indexes_list)
        q_m = self._mean(q)
        q_v = self._variance_activation(self._variance(q)) + self.variance_eps
        z_dist = self._distribution(q_m, torch.sqrt(q_v))
        z = self.z_transformation(z_dist.rsample())

        if self.return_distribution:
            return z_dist, z
        return q_m, q_v, z
