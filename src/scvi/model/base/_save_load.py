from __future__ import annotations

import inspect
import logging
import os
import warnings
from typing import TYPE_CHECKING

import anndata
import mudata
import numpy as np
import pandas as pd
import torch
from anndata import AnnData, read_h5ad

from scvi import settings
from scvi.data._constants import _SETUP_METHOD_NAME
from scvi.data._download import _download
from scvi.model.base._constants import SAVE_KEYS

if TYPE_CHECKING:
    from typing import Literal

    import numpy.typing as npt

    from scvi._types import AnnOrMuData

logger = logging.getLogger(__name__)


def _load_legacy_saved_files(
    dir_path: str,
    file_name_prefix: str,
    load_adata: bool,
) -> tuple[dict, np.ndarray, dict, AnnData | None]:
    model_path = os.path.join(dir_path, f"{file_name_prefix}{SAVE_KEYS.LEGACY_MODEL_FNAME}")
    var_names_path = os.path.join(
        dir_path, f"{file_name_prefix}{SAVE_KEYS.LEGACY_VAR_NAMES_FNAME}"
    )
    setup_dict_path = os.path.join(
        dir_path, f"{file_name_prefix}{SAVE_KEYS.LEGACY_SETUP_DICT_FNAME}"
    )

    model_state_dict = torch.load(model_path, map_location="cpu", weights_only=False)

    var_names = np.genfromtxt(var_names_path, delimiter=",", dtype=str)

    with open(setup_dict_path, "rb") as handle:
        attr_dict = pd.read_pickle(handle)

    if load_adata:
        adata_path = os.path.join(dir_path, f"{file_name_prefix}{SAVE_KEYS.ADATA_FNAME}")
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
    prefix: str | None = None,
    map_location: Literal["cpu", "cuda"] | None = None,
    backup_url: str | None = None,
) -> tuple[dict, np.ndarray, dict, AnnData]:
    """Helper to load saved files."""
    file_name_prefix = prefix or ""

    model_file_name = f"{file_name_prefix}{SAVE_KEYS.MODEL_FNAME}"
    model_path = os.path.join(dir_path, model_file_name)
    try:
        _download(backup_url, dir_path, model_file_name)
        model = torch.load(model_path, map_location=map_location, weights_only=False)
    except FileNotFoundError as exc:
        raise ValueError(
            f"Failed to load model file at {model_path}. "
            "If attempting to load a saved model from <v0.15.0, please use the util function "
            "`convert_legacy_save` to convert to an updated format."
        ) from exc

    model_state_dict = model.get(SAVE_KEYS.MODEL_STATE_DICT_KEY)
    var_names = model.get(SAVE_KEYS.VAR_NAMES_KEY)
    attr_dict = model.get(SAVE_KEYS.ATTR_DICT_KEY)

    if load_adata:
        is_mudata = attr_dict["registry_"].get(_SETUP_METHOD_NAME) == "setup_mudata"
        file_suffix = SAVE_KEYS.ADATA_FNAME if is_mudata is False else SAVE_KEYS.MDATA_FNAME
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


def _initialize_model(cls, adata, registry, attr_dict, datamodule):
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

    if not adata:
        adata = None

    if datamodule:
        non_kwargs["datamodule"] = datamodule

    if "registry" in inspect.signature(cls).parameters:
        model = cls(adata, registry=registry, **non_kwargs, **kwargs)
    else:
        model = cls(adata, **non_kwargs, **kwargs)
    for attr, val in attr_dict.items():
        setattr(model, attr, val)

    return model


def _get_var_names(
    adata: AnnOrMuData,
    legacy_mudata_format: bool = False,
) -> npt.NDArray | dict[str, npt.NDArray]:
    """Get variable names from an :class:`~anndata.AnnData` or :class:`~mudata.MuData` object.

    Parameters
    ----------
    adata
        :class:`~anndata.AnnData` or :class:`~mudata.MuData` object.
    legacy_mudata_format
        If ``True``, returns variable names for :class:`~mudata.MuData` objects in the legacy
        format, _i.e._, a flat array with variable names across all modalities concatenated
        together. Otherwise, returns a dictionary with keys corresponding to the modality names
        and values corresponding to the variable names for each modality.

    Returns
    -------
    - An array of variable names if ``adata`` is an :class:`~anndata.AnnData` object
    - An array of concatenated modality variable names if ``adata`` is a :class:`~mudata.MuData`
        and ``legacy_mudata_format`` is ``True``
    - A dictionary with keys corresponding to the modality names and values corresponding to the
        variable names for each modality if ``adata`` is a :class:`~mudata.MuData` and
        ``legacy_mudata_format`` is ``False``
    """
    if isinstance(adata, AnnData) or legacy_mudata_format:
        return adata.var_names.astype(str).to_numpy()
    elif isinstance(adata, mudata.MuData):
        return {
            mod_key: adata.mod[mod_key].var_names.astype(str).to_numpy()
            for mod_key in adata.mod.keys()
        }
    else:
        raise TypeError("`adata` must be an AnnData or MuData object.")


def _validate_var_names(
    adata: AnnOrMuData | None,
    source_var_names: npt.NDArray | dict[str, npt.NDArray],
    load_var_names: npt.NDArray | dict[str, npt.NDArray] | None = None,
) -> None:
    """Validate that source and loaded variable names match.

    Parameters
    ----------
    adata
        :class:`~anndata.AnnData` or :class:`~mudata.MuData` object.
    source_var_names
        Variable names from a saved model file corresponding to the variable names used during
        training.
    load_var_names
        Variable names from the loaded registry.
    """
    from numpy import array_equal

    source_per_mod_var_names = isinstance(source_var_names, dict)

    if load_var_names is None:
        is_anndata = isinstance(adata, AnnData)
        load_var_names = _get_var_names(
            adata,
            legacy_mudata_format=(not is_anndata and not source_per_mod_var_names),
        )

    if source_per_mod_var_names:
        valid_load_var_names = all(
            array_equal(source_var_names.get(mod_key), load_var_names.get(mod_key))
            for mod_key in source_var_names
        )
    else:
        valid_load_var_names = array_equal(source_var_names, load_var_names)

    if not valid_load_var_names:
        warnings.warn(
            "`var_names` for the loaded `model` does not match those used to "
            "train the model. For valid results, the former should match the latter.",
            UserWarning,
            stacklevel=settings.warnings_stacklevel,
        )
