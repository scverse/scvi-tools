from typing import List, Literal, Optional

import torch
from chex import dataclass
from torch import nn

from ._utils import _one_hot_torch

BATCH_NORM_KWARGS = {
    "momentum": 0.01,
    "epsilon": 0.001,
}
LAYER_NORM_KWARGS = {
    "elementwise_affine": False,
}


@dataclass
class Module(nn.Module):
    """Class inheriting from :class:`~torch.nn.Module` and :class:`~chex.dataclass`."""

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


class TorchFCLayers(Module):
    """Fully-connected layers with a PyTorch backend."""

    n_input: int
    layers_dim: List[int]
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

    def setup(self):
        """Setup the module."""
        self._norm = None
        if self.norm == "batch":
            self._norm = BatchNorm
            self._norm_kwargs = BATCH_NORM_KWARGS.copy()
            self._norm_kwargs.update(self.norm_kwargs)
        elif self.norm == "layer":
            self._norm = nn.LayerNorm
            self._norm_kwargs = LAYER_NORM_KWARGS.copy()
            self._norm_kwargs.update(self.norm_kwargs)
        elif self.norm == "group":
            self._norm = nn.GroupNorm

        self._activation = None
        if self.activation:
            self._activation = getattr(nn, self.activation)

        _cat_dim = sum(self.n_classes_list)
        _emb_dim = self.embedding_dim or _cat_dim
        _layers_dim = [self.n_input] + self.layers_dim
        self._module = nn.ModuleList()
        for i, (n_in, n_out) in enumerate(zip(_layers_dim[:-1], _layers_dim[1:])):
            inject = self._inject_into_layer(i)
            layer = []
            layer.append(
                Linear(
                    n_in + _emb_dim if inject else n_in,
                    n_out,
                    bias=self.bias,
                    inject_embedding=inject,
                )
            )
            if self._norm:
                layer.append(self._norm(n_out, **self._norm_kwargs))
            if self._activation:
                layer.append(self._activation(**self.activation_kwargs))
            if self.dropout_rate:
                layer.append(nn.Dropout(p=self.dropout_rate))
            self._module.append(nn.Sequential(*layer))

        self._embedding = None
        if self.embedding_dim:
            # when embedding_dim = None, use one-hot encoding instead
            self._embedding = nn.Embedding(_cat_dim, self.embedding_dim)

    def _inject_into_layer(self, layer_num: int) -> bool:
        """Whether to inject covariates into a given layer."""
        return layer_num == 0 or self.inject_covariates

    def set_online_update_hooks(self, hook_first_layer: bool = True):
        """Set online update hooks for transfer learning."""
        self._hooks = []

        def _hook_fn_weight(grad):
            cat_dim = sum(self.n_classes_list)
            new_grad = torch.zeros_like(grad)
            if cat_dim > 0:
                new_grad[:, -cat_dim:] = grad[:, -cat_dim:]
            return new_grad

        def _hook_fn_zero_out(grad):
            return grad * 0

        for i, block in enumerate(self._module):
            for layer in block:
                if (i == 0 and not hook_first_layer) or (not isinstance(layer, Linear)):
                    continue
                hook = _hook_fn_weight if layer._inject_embedding else _hook_fn_zero_out
                w = layer.weight.register_hook(hook)
                b = layer.bias.register_hook(_hook_fn_zero_out)
                self._hooks += [w, b]

    def _get_covariates(self, indexes_list: List[torch.Tensor]) -> torch.Tensor:
        """
        Converts a list of category membership indexes into an embedding.

        Embedding is either a one-hot encoded :class:`~torch.Tensor` or an embedding
        :class:`~torch.Tensor`.

        Parameters
        ----------
        indexes_list
            List of :class:`~torch.Tensor` of shape `(n_samples, 1)` or
            `(n_samples, n_classes)` containing the categorical memberships of each
            sample for a given covariate, where the order of covariates is the same as
            in `self.n_classes_list`.
        """
        one_hot_list = []
        for n_classes, indexes in zip(self.n_classes_list, indexes_list):
            if n_classes > 0 and indexes is None:
                raise ValueError(
                    "Indexes must be provided for covariates with `n_classes` > 0."
                )
            if n_classes <= 1:  # n_cats = 0,1 adds no information
                continue
            if indexes.size(-1) == n_classes:  # already one-hot encoded
                one_hot = indexes
            else:
                one_hot = _one_hot_torch(indexes, n_classes=n_classes)
            one_hot_list.append(one_hot)
        covariates = torch.cat(one_hot_list, dim=-1)  # (n_samples, n_total_classes)

        if self._embedding is not None:
            indexes = torch.nonzero(one_hot)[:, -1]  # (n_samples, 1)
            covariates = self._embedding(indexes)  # (n_samples, embedding_dim)

        return covariates

    def forward(
        self,
        x: torch.Tensor,
        indexes_list: Optional[List[torch.Tensor]] = None,
    ) -> torch.Tensor:
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
            sample for a given covariate, where the order of covariates is the same as
            in `self.n_classes_list`.

        Returns
        -------
        out
            Output :class:`~torch.Tensor` of shape `(n_samples, n_output)` or
            `(n_layers, n_samples, n_output)`.
        """
        if len(self.n_cat_list) > len(indexes_list):
            raise ValueError("Not enough covariates provided in `indexes_list`.")
        covariates = None
        if indexes_list is not None:
            covariates = self._get_covariates(indexes_list)

        out = x
        for block in self._module:
            x_ = out
            for layer in block:
                if isinstance(layer, Linear):
                    # only inject covariates into linear layers
                    out = layer(out, embedding=covariates)
                else:
                    out = layer(out)
            if self.residual:
                out = out + x_
        return out
