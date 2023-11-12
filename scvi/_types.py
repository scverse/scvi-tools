from typing import Literal, Union

import anndata
import mudata
import torch

from scvi._packageproxy import jnp

Number = Union[int, float]
AnnOrMuData = Union[anndata.AnnData, mudata.MuData]
try:
    Tensor = Union[torch.Tensor, jnp.ndarray]
except ImportError:
    Tensor = torch.Tensor
LossRecord = Union[dict[str, Tensor], Tensor]
# TODO(adamgayoso): Add constants for minified data types.
MinifiedDataType = Literal["latent_posterior_parameters"]
