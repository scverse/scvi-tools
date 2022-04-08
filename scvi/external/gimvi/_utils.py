import os
import pickle
from typing import Optional, Tuple

import numpy as np
import torch
from anndata import AnnData, read

from scvi._compat import Literal
from scvi.model._utils import _download_if_missing


def _load_legacy_saved_gimvi_files(
    dir_path: str,
    file_name_prefix: str,
    load_seq_adata: bool,
    load_spatial_adata: bool,
    backup_url: Optional[str] = None,
) -> Tuple[dict, np.ndarray, np.ndarray, dict, Optional[AnnData], Optional[AnnData]]:
    model_path = os.path.join(dir_path, f"{file_name_prefix}model_params.pt")
    setup_dict_path = os.path.join(dir_path, f"{file_name_prefix}attr.pkl")
    seq_var_names_path = os.path.join(dir_path, f"{file_name_prefix}var_names_seq.csv")
    spatial_var_names_path = os.path.join(
        dir_path, f"{file_name_prefix}var_names_spatial.csv"
    )

    _download_if_missing(model_path, backup_url)
    model_state_dict = torch.load(model_path)

    _download_if_missing(seq_var_names_path, backup_url)
    seq_var_names = np.genfromtxt(seq_var_names_path, delimiter=",", dtype=str)
    _download_if_missing(spatial_var_names_path, backup_url)
    spatial_var_names = np.genfromtxt(spatial_var_names_path, delimiter=",", dtype=str)

    _download_if_missing(setup_dict_path, backup_url)
    with open(setup_dict_path, "rb") as handle:
        attr_dict = pickle.load(handle)

    adata_seq, adata_spatial = None, None
    if load_seq_adata:
        seq_data_path = os.path.join(dir_path, f"{file_name_prefix}adata_seq.h5ad")
        _download_if_missing(seq_data_path, backup_url)
        if os.path.exists(seq_data_path):
            adata_seq = read(seq_data_path)
        elif not os.path.exists(seq_data_path):
            raise ValueError(
                "Save path contains no saved anndata and no adata was passed."
            )
    if load_spatial_adata:
        spatial_data_path = os.path.join(
            dir_path, f"{file_name_prefix}adata_spatial.h5ad"
        )
        _download_if_missing(spatial_data_path, backup_url)
        if os.path.exists(spatial_data_path):
            adata_spatial = read(spatial_data_path)
        elif not os.path.exists(spatial_data_path):
            raise ValueError(
                "Save path contains no saved anndata and no adata was passed."
            )

    return (
        model_state_dict,
        seq_var_names,
        spatial_var_names,
        attr_dict,
        adata_seq,
        adata_spatial,
    )


def _load_saved_gimvi_files(
    dir_path: str,
    load_seq_adata: bool,
    load_spatial_adata: bool,
    prefix: Optional[str] = None,
    map_location: Optional[Literal["cpu", "cuda"]] = None,
    backup_url: Optional[str] = None,
) -> Tuple[
    dict, dict, np.ndarray, np.ndarray, dict, Optional[AnnData], Optional[AnnData]
]:
    file_name_prefix = prefix or ""

    model_path = os.path.join(dir_path, f"{file_name_prefix}model.pt")
    try:
        _download_if_missing(model_path, backup_url)
        model = torch.load(model_path, map_location=map_location)
    except Exception as exc:
        raise ValueError(
            f"Failed to load model file at {model_path}. "
            "If attempting to load a saved model from <v0.15.0, please use the util function "
            "`convert_legacy_save` to convert to an updated format."
        ) from exc

    model_state_dict = model["model_state_dict"]
    seq_var_names = model["seq_var_names"]
    spatial_var_names = model["spatial_var_names"]
    attr_dict = model["attr_dict"]

    adata_seq, adata_spatial = None, None
    if load_seq_adata:
        seq_data_path = os.path.join(dir_path, f"{file_name_prefix}adata_seq.h5ad")
        _download_if_missing(seq_data_path, backup_url)
        if os.path.exists(seq_data_path):
            adata_seq = read(seq_data_path)
        elif not os.path.exists(seq_data_path):
            raise ValueError(
                "Save path contains no saved anndata and no adata was passed."
            )
    if load_spatial_adata:
        spatial_data_path = os.path.join(
            dir_path, f"{file_name_prefix}adata_spatial.h5ad"
        )
        _download_if_missing(spatial_data_path, backup_url)
        if os.path.exists(spatial_data_path):
            adata_spatial = read(spatial_data_path)
        elif not os.path.exists(spatial_data_path):
            raise ValueError(
                "Save path contains no saved anndata and no adata was passed."
            )

    return (
        attr_dict,
        seq_var_names,
        spatial_var_names,
        model_state_dict,
        adata_seq,
        adata_spatial,
    )
