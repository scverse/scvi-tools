from __future__ import annotations

from typing import Literal

import anndata
import mudata
import torch

from scvi.utils import is_package_installed

Number = int | float
AnnOrMuData = anndata.AnnData | mudata.MuData
if is_package_installed("jax"):
    import jax.numpy as jnp

    Tensor = torch.Tensor | jnp.ndarray
else:
    Tensor = torch.Tensor
LossRecord = dict[str, Tensor] | Tensor
MinifiedDataType = Literal["latent_posterior_parameters"]
