import logging
import os
import warnings
from typing import Literal, Optional

import anndata
import mudata
import numpy as np
import pandas as pd
import torch
from anndata import AnnData, read_h5ad

from scvi import settings
from scvi.data._constants import _SETUP_METHOD_NAME
from scvi.data._download import _download

logger = logging.getLogger(__name__)


def _load_legacy_saved_files(
    dir_path: str,
    file_name_prefix: str,
    load_adata: bool,
) -> tuple[dict, np.ndarray, dict, Optional[AnnData]]:
    model_path = os.path.join(dir_path, f"{file_name_prefix}model_params.pt")
    var_names_path = os.path.join(dir_path, f"{file_name_prefix}var_names.csv")
    setup_dict_path = os.path.join(dir_path, f"{file_name_prefix}attr.pkl")

    model_state_dict = torch.load(model_path, map_location="cpu")

    var_names = np.genfromtxt(var_names_path, delimiter=",", dtype=str)

    with open(setup_dict_path, "rb") as handle:
        attr_dict = pd.read_pickle(handle)

    if load_adata:
        adata_path = os.path.join(dir_path, f"{file_name_prefix}adata.h5ad")
        if os.path.exists(adata_path):
            adata = read_h5ad(adata_path)
        elif not os.path.exists(adata_path):
            raise ValueError("Save path contains no saved anndata and no adata was passed.")
    else:
        adata = None

    return model_state_dict, var_names, attr_dict, adata


def _load_saved_files(
    dir_path: str,
    load_adata: bool,
    prefix: Optional[str] = None,
    map_location: Optional[Literal["cpu", "cuda"]] = None,
    backup_url: Optional[str] = None,
) -> tuple[dict, np.ndarray, dict, AnnData]:
    """Helper to load saved files."""
    file_name_prefix = prefix or ""

    model_file_name = f"{file_name_prefix}model.pt"
    model_path = os.path.join(dir_path, model_file_name)
    try:
        _download(backup_url, dir_path, model_file_name)
        model = torch.load(model_path, map_location=map_location)
    except FileNotFoundError as exc:
        raise ValueError(
            f"Failed to load model file at {model_path}. "
            "If attempting to load a saved model from <v0.15.0, please use the util function "
            "`convert_legacy_save` to convert to an updated format."
        ) from exc

    model_state_dict = model["model_state_dict"]
    var_names = model["var_names"]
    attr_dict = model["attr_dict"]

    if load_adata:
        is_mudata = attr_dict["registry_"].get(_SETUP_METHOD_NAME) == "setup_mudata"
        file_suffix = "adata.h5ad" if is_mudata is False else "mdata.h5mu"
        adata_path = os.path.join(dir_path, f"{file_name_prefix}{file_suffix}")
        if os.path.exists(adata_path):
            if is_mudata:
                adata = mudata.read(adata_path)
            else:
                adata = anndata.read_h5ad(adata_path)
        else:
            raise ValueError("Save path contains no saved anndata and no adata was passed.")
    else:
        adata = None

    return attr_dict, var_names, model_state_dict, adata


def _initialize_model(cls, adata, attr_dict):
    """Helper to initialize a model."""
    if "init_params_" not in attr_dict.keys():
        raise ValueError(
            "No init_params_ were saved by the model. Check out the "
            "developers guide if creating custom models."
        )
    # get the parameters for the class init signature
    init_params = attr_dict.pop("init_params_")

    # new saving and loading, enable backwards compatibility
    if "non_kwargs" in init_params.keys():
        # grab all the parameters except for kwargs (is a dict)
        non_kwargs = init_params["non_kwargs"]
        kwargs = init_params["kwargs"]

        # expand out kwargs
        kwargs = {k: v for (i, j) in kwargs.items() for (k, v) in j.items()}
    else:
        # grab all the parameters except for kwargs (is a dict)
        non_kwargs = {k: v for k, v in init_params.items() if not isinstance(v, dict)}
        kwargs = {k: v for k, v in init_params.items() if isinstance(v, dict)}
        kwargs = {k: v for (i, j) in kwargs.items() for (k, v) in j.items()}
        non_kwargs.pop("use_cuda")

    # backwards compat for scANVI
    if "unlabeled_category" in non_kwargs.keys():
        non_kwargs.pop("unlabeled_category")
    if "pretrained_model" in non_kwargs.keys():
        non_kwargs.pop("pretrained_model")

    model = cls(adata, **non_kwargs, **kwargs)
    for attr, val in attr_dict.items():
        setattr(model, attr, val)

    return model


def _validate_var_names(adata, source_var_names):
    user_var_names = adata.var_names.astype(str)
    if not np.array_equal(source_var_names, user_var_names):
        warnings.warn(
            "var_names for adata passed in does not match var_names of adata used to "
            "train the model. For valid results, the vars need to be the same and in "
            "the same order as the adata used to train the model.",
            UserWarning,
            stacklevel=settings.warnings_stacklevel,
        )
