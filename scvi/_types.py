from inspect import isfunction
from typing import Any, Literal, Union

import anndata
import jax.numpy as jnp
import mudata
import torch

from scvi._decorators import classproperty

Number = Union[int, float]
AnnOrMuData = Union[anndata.AnnData, mudata.MuData]
Tensor = Union[torch.Tensor, jnp.ndarray]
LossRecord = Union[dict[str, Tensor], Tensor]
# TODO(adamgayoso): Add constants for minified data types.
MinifiedDataType = Literal["latent_posterior_parameters"]


class TunableMeta(type):
    """Metaclass for Tunable class."""

    def __getitem__(cls, values):
        if not isinstance(values, tuple):
            values = (values,)
        return type("Tunable_", (Tunable,), {"__args__": values})


class Tunable(metaclass=TunableMeta):
    """Typing class for tagging keyword arguments as tunable."""


class TunableMixin:
    """Mixin class for exposing tunable attributes."""

    @classproperty
    def _tunables(cls) -> list[Any]:
        """Returns the tunable attributes of the model class."""
        _tunables = []
        for attr_key in dir(cls):
            if attr_key == "_tunables":
                # Don't recurse
                continue
            attr = getattr(cls, attr_key)
            if hasattr(attr, "_tunables") or isfunction(attr):
                _tunables.append(attr)
        return _tunables
