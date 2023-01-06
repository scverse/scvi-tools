from typing import List, Literal, Optional

import torch
from chex import dataclass
from torch import nn
from torch.nn import functional as F

from ._utils import one_hot

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


class TorchFCLayers(Module):
    """Fully-connected layers with a PyTorch backend."""

    layers_dim: List[int]
    bias: bool
    dropout_rate: float
    norm: Literal["batch", "layer", "group", None]
    norm_kwargs: dict
    activation: str
    activation_kwargs: dict
    residual: bool
    n_cat_list: List[int]
    inject_covariates: bool
    training: Optional[bool]

    def setup(self):
        """Setup the module."""
        self._norm = None
        if self.norm == "batch":
            self._norm = nn.BatchNorm1d
            self._norm_kwargs = BATCH_NORM_KWARGS.copy()
            self._norm_kwargs.update(self.norm_kwargs)
        elif self.norm == "layer":
            self._norm = nn.LayerNorm
            self._norm_kwargs = LAYER_NORM_KWARGS.copy()
            self._norm_kwargs.update(self.norm_kwargs)
        elif self.norm == "group":
            self._norm = nn.GroupNorm

        self._activation = getattr(F, self.activation) if self.activation else None

        cat_dim = sum(self.n_cat_list)
        self._module = nn.ModuleList()
        for i, (n_in, n_out) in enumerate(
            zip(self.layers_dim[:-1], self.layers_dim[1:])
        ):
            if self._inject_into_layer(i):
                n_in += cat_dim
            layer = [nn.Linear(n_in, n_out, bias=self.bias)]
            if self._norm:
                layer.append(self._norm(n_out, **self._norm_kwargs))
            # activation and dropout not included - use function versions
            self._module.append(nn.Sequential(*layer))

    def _inject_into_layer(self, layer_num: int) -> bool:
        """Whether to inject covariates into a given layer."""
        return layer_num == 0 or self.inject_covariates

    def forward(
        self,
        x: torch.Tensor,
        cat_list: Optional[List[torch.Tensor]] = None,
    ) -> torch.Tensor:
        """
        Forward pass through the module.

        Parameters
        ----------
        x
            Input tensor of shape (batch_size, n_in).
        cat_list
            List of categorical membership(s) for each sample.

        Returns
        -------
        out
            Output :class:`~torch.Tensor` of shape (batch_size, n_out).
        """
        one_hot_cat_list = []

        if len(self.n_cat_list) > len(cat_list):
            raise ValueError("Categorical arguments provided don't match init params.")
        for n_cat, cat in zip(self.n_cat_list, cat_list):
            if n_cat and cat is None:
                raise ValueError("`cat` not provided while `n_cat` is not 0.")
            if n_cat == 1:
                continue
            if cat.size(1) != n_cat:
                one_hot_cat = one_hot(cat, n_cat)
            else:
                one_hot_cat = cat
            one_hot_cat_list += [one_hot_cat]
