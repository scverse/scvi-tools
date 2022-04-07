import os
import pickle
from typing import Optional, Tuple

import numpy as np
import torch
from anndata import AnnData, read
from sklearn.utils import deprecated

from scvi._compat import Literal


def _should_use_legacy_saved_gimvi_files(dir_path: str, file_name_prefix: str) -> bool:
    new_model_path = os.path.join(dir_path, f"{file_name_prefix}model.pt")
    model_path = os.path.join(dir_path, f"{file_name_prefix}model_params.pt")
    seq_var_names_path = os.path.join(dir_path, f"{file_name_prefix}var_names_seq.csv")
    spatial_var_names_path = os.path.join(
        dir_path, f"{file_name_prefix}var_names_spatial.csv"
    )
    attr_dict_path = os.path.join(dir_path, f"{file_name_prefix}attr.pkl")
    return (
        not os.path.exists(new_model_path)
        and os.path.exists(model_path)
        and os.path.exists(seq_var_names_path)
        and os.path.exists(spatial_var_names_path)
        and os.path.exists(attr_dict_path)
    )


@deprecated(
    extra="Please update your saved models to use the latest version. The legacy save and load scheme will be removed in version 0.16.0."
)
def _load_legacy_saved_gimvi_files(
    dir_path: str,
    file_name_prefix: str,
    map_location: Optional[Literal["cpu", "cuda"]],
) -> Tuple[dict, np.ndarray, np.ndarray, dict]:
    model_path = os.path.join(dir_path, f"{file_name_prefix}model_params.pt")
    setup_dict_path = os.path.join(dir_path, f"{file_name_prefix}attr.pkl")
    seq_var_names_path = os.path.join(dir_path, f"{file_name_prefix}var_names_seq.csv")
    spatial_var_names_path = os.path.join(
        dir_path, f"{file_name_prefix}var_names_spatial.csv"
    )

    model_state_dict = torch.load(model_path, map_location=map_location)

    seq_var_names = np.genfromtxt(seq_var_names_path, delimiter=",", dtype=str)
    spatial_var_names = np.genfromtxt(spatial_var_names_path, delimiter=",", dtype=str)

    with open(setup_dict_path, "rb") as handle:
        attr_dict = pickle.load(handle)

    return model_state_dict, seq_var_names, spatial_var_names, attr_dict


def _load_saved_gimvi_files(
    dir_path: str,
    load_seq_adata: bool,
    load_spatial_adata: bool,
    prefix: Optional[str] = None,
    map_location: Optional[Literal["cpu", "cuda"]] = None,
) -> Tuple[dict, dict, np.ndarray, np.ndarray, dict, AnnData, AnnData]:
    file_name_prefix = prefix or ""
    seq_data_path = os.path.join(dir_path, f"{file_name_prefix}adata_seq.h5ad")
    spatial_data_path = os.path.join(dir_path, f"{file_name_prefix}adata_spatial.h5ad")

    adata_seq, adata_spatial = None, None
    if load_seq_adata and os.path.exists(seq_data_path):
        adata_seq = read(seq_data_path)
    elif load_seq_adata and not os.path.exists(seq_data_path):
        raise ValueError("Save path contains no saved anndata and no adata was passed.")
    if load_spatial_adata and os.path.exists(spatial_data_path):
        adata_spatial = read(spatial_data_path)
    elif load_spatial_adata and not os.path.exists(spatial_data_path):
        raise ValueError("Save path contains no saved anndata and no adata was passed.")

    use_legacy = _should_use_legacy_saved_gimvi_files(dir_path, file_name_prefix)

    # TODO(jhong): Remove once legacy load is deprecated.
    if use_legacy:
        (
            model_state_dict,
            seq_var_names,
            spatial_var_names,
            attr_dict,
        ) = _load_legacy_saved_gimvi_files(dir_path, file_name_prefix, map_location)
    else:
        model_path = os.path.join(dir_path, f"{file_name_prefix}model.pt")

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
