from typing import Dict, Union

import anndata
import jax.numpy as jnp
import mudata
import torch

from scvi._compat import Literal

Number = Union[int, float]
AnnOrMuData = Union[anndata.AnnData, mudata.MuData]
Tensor = Union[torch.Tensor, jnp.ndarray]
LossRecord = Union[Dict[str, Tensor], Tensor]
LatentDataType = Union[Literal["posterior_parameters"], None]
