import os
from typing import Optional, Tuple

import numpy as np
import torch
from anndata import AnnData, read

from scvi._compat import Literal
from scvi.model._utils import _download_if_missing


def _load_saved_gimvi_files(
    dir_path: str,
    load_seq_adata: bool,
    load_spatial_adata: bool,
    prefix: Optional[str] = None,
    map_location: Optional[Literal["cpu", "cuda"]] = None,
    backup_url: Optional[str] = None,
) -> Tuple[dict, dict, np.ndarray, np.ndarray, dict, AnnData, AnnData]:
    file_name_prefix = prefix or ""
    seq_data_path = os.path.join(dir_path, f"{file_name_prefix}adata_seq.h5ad")
    spatial_data_path = os.path.join(dir_path, f"{file_name_prefix}adata_spatial.h5ad")

    adata_seq, adata_spatial = None, None
    if load_seq_adata:
        _download_if_missing(seq_data_path, backup_url)
        if os.path.exists(seq_data_path):
            adata_seq = read(seq_data_path)
        elif not os.path.exists(seq_data_path):
            raise ValueError(
                "Save path contains no saved anndata and no adata was passed."
            )
    if load_spatial_adata:
        _download_if_missing(spatial_data_path, backup_url)
        if os.path.exists(spatial_data_path):
            adata_spatial = read(spatial_data_path)
        elif not os.path.exists(spatial_data_path):
            raise ValueError(
                "Save path contains no saved anndata and no adata was passed."
            )

    model_path = os.path.join(dir_path, f"{file_name_prefix}model.pt")
    _download_if_missing(model_path, backup_url)

    model = torch.load(model_path, map_location=map_location)
    model_state_dict = model["model_state_dict"]
    seq_var_names = model["seq_var_names"]
    spatial_var_names = model["spatial_var_names"]
    attr_dict = model["attr_dict"]

    return (
        attr_dict,
        seq_var_names,
        spatial_var_names,
        model_state_dict,
        adata_seq,
        adata_spatial,
    )
