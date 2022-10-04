from typing import Dict, Union

import anndata
import jax.numpy as jnp
import mudata
import torch

Number = Union[int, float]
AnnOrMuData = Union[anndata.AnnData, mudata.MuData]
Tensor = Union[torch.Tensor, jnp.ndarray]
LossRecord = Union[Dict[str, Tensor], Tensor]
LatentDataType = str
