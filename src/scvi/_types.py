from __future__ import annotations

from typing import Literal

import anndata
import mudata
import torch

Number = int | float
AnnOrMuData = anndata.AnnData | mudata.MuData
Tensor = torch.Tensor
LossRecord = dict[str, Tensor] | Tensor
MinifiedDataType = Literal["latent_posterior_parameters"]
