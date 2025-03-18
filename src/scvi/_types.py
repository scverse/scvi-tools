from __future__ import annotations

from typing import Literal

import anndata
import jax.numpy as jnp
import mudata
import torch

Number = int | float
AnnOrMuData = anndata.AnnData | mudata.MuData
Tensor = torch.Tensor | jnp.ndarray
LossRecord = dict[str, Tensor] | Tensor
# TODO(adamgayoso): Add constants for minified data types.
MinifiedDataType = Literal["latent_posterior_parameters"]
