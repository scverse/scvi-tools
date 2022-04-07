from typing import Dict, Type, Union

import jax.numpy as jnp
import torch

from scvi.data.fields import BaseAnnDataField

Number = Union[int, float]
AnnDataField = Type[BaseAnnDataField]
Tensor = Union[torch.Tensor, jnp.ndarray]
LossRecord = Union[Dict[str, Tensor], Tensor]
