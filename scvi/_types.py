from typing import Dict, Type, Union

import anndata
import jax.numpy as jnp
import mudata
import torch

from scvi.data.fields import BaseAnnDataField

Number = Union[int, float]
AnnOrMuData = Union[anndata.AnnData, mudata.MuData]
AnnDataField = Type[BaseAnnDataField]
Tensor = Union[torch.Tensor, jnp.ndarray]
LossRecord = Union[Dict[str, Tensor], Tensor]
