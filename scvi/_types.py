from typing import Literal, Union

import anndata
import jax.numpy as jnp
import mudata
import torch

Number = Union[int, float]
AnnOrMuData = Union[anndata.AnnData, mudata.MuData]
Tensor = Union[torch.Tensor, jnp.ndarray]
LossRecord = Union[dict[str, Tensor], Tensor]
# TODO(adamgayoso): Add constants for minified data types.
MinifiedDataType = Literal["latent_posterior_parameters"]
