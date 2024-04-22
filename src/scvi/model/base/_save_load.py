from __future__ import annotations

import logging
import os
import warnings
from typing import Any, Literal

import mudata
import numpy as np
import numpy.typing as npt
import pandas as pd
import torch
from anndata import AnnData, read_h5ad
from torch import Tensor

from scvi import settings
from scvi._types import AnnOrMuData
from scvi.model.base._constants import SAVE_KEYS

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

    model_state_dict = torch.load(model_path, map_location="cpu")

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
) -> tuple[
    dict[str, Any], npt.NDArray | dict[str, npt.NDArray], dict[str, Tensor], AnnOrMuData | None
]:
    """Load saved model files and data from a directory.

    Parameters
    ----------
    dir_path
        Directory path containing the saved model files.
    load_adata
        Whether to attempt to load a saved ``h5ad`` or ``h5mu`` file. If ``True`` and no such file
        exists, a ``FileNotFoundError`` is raised.

    Returns
    -------
    A tuple containing the following elements:

    * Dictionary of model attributes.
    * ``var_names`` from the :class:`~anndata.AnnData` or :class:`~mudata.MuData` object used to
        train the model.
    * Model state dictionary mapping parameter names to :class:`~torch.Tensor`s.
    * Loaded :class:`~anndata.AnnData` or :class:`~mudata.MuData` object if ``load_adata`` is
        ``True``, else ``None``.

    Raises
    ------
    ``ValueError``
        If the model file cannot be loaded.
    ``FileNotFoundError``
        If ``load_adata`` is ``True`` and no saved ``h5ad`` or ``h5mu`` file exists in
        ``dir_path``.
    """
    from anndata import read_h5ad
    from mudata import read_h5mu

    from scvi.data._constants import _SETUP_METHOD_NAME
    from scvi.data._download import _download

    prefix = prefix or ""

    model_file_name = f"{prefix}{SAVE_KEYS.MODEL_FNAME}"
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

    model_state_dict = model.get(SAVE_KEYS.MODEL_STATE_DICT_KEY)
    var_names = model.get(SAVE_KEYS.VAR_NAMES_KEY)
    attr_dict = model.get(SAVE_KEYS.ATTR_DICT_KEY)

    is_anndata = attr_dict.get("registry_").get(_SETUP_METHOD_NAME) == "setup_anndata"
    file_suffix = SAVE_KEYS.ADATA_FNAME if is_anndata else SAVE_KEYS.MDATA_FNAME
    read_method = read_h5ad if is_anndata else read_h5mu
    data_path = os.path.join(dir_path, f"{prefix}{file_suffix}")

    if load_adata and not os.path.exists(data_path):
        raise FileNotFoundError(
            "Save path contains no saved h5ad or h5mu file and no adata was passed."
        )
    elif load_adata:
        adata = read_method(data_path)
    else:
        adata = None

    return attr_dict, var_names, model_state_dict, adata


def _initialize_model(cls, adata: AnnOrMuData, attr_dict: dict[str, Any]):
    """Helper to initialize a model."""
    if "init_params_" not in attr_dict:
        raise ValueError(
            "No init_params_ were saved by the model. Check out the "
            "developers guide if creating custom models."
        )
    # get the parameters for the class init signature
    init_params = attr_dict.pop("init_params_")

    # new saving and loading, enable backwards compatibility
    if "non_kwargs" in init_params:
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
    if "unlabeled_category" in non_kwargs:
        non_kwargs.pop("unlabeled_category")
    if "pretrained_model" in non_kwargs:
        non_kwargs.pop("pretrained_model")

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
        Whether to return variable names for :class:`~mudata.MuData` objects in the legacy format,
        where variable names across modalities are concatenated into a single array.

    Returns
    -------
    An array of variable names if ``adata`` is an :class:`~anndata.AnnData` object, or a
    dictionary of variable names if ``adata`` is a :class:`~mudata.MuData` object, where keys
    correspond to modalities and values are arrays of variable names.
    """
    if isinstance(adata, AnnData) or legacy_mudata_format:
        return adata.var_names.astype(str).to_numpy()
    elif isinstance(adata, mudata.MuData):
        return {
            mod_key: adata.mod[mod_key].var_names.astype(str).to_numpy() for mod_key in adata.mod
        }
    else:
        raise TypeError("`adata` must be an AnnData or MuData object.")


def _validate_var_names(
    adata: AnnData, source_var_names: npt.NDArray | dict[str, npt.NDArray]
) -> None:
    """Validate that source and loaded ``var_names`` are the same.

    Parameters
    ----------
    adata
        :class:`~anndata.AnnData` or :class:`~mudata.MuData` object.
    source_var_names
        ``var_names`` from a saved model file corresponding to the variable names used during
        model training.

    Raises
    ------
    UserWarning
        If ``var_names`` for the loaded ``adata`` do not match those of the ``adata`` used to train
        the model.
    """
    from numpy import array_equal

    is_anndata = isinstance(adata, AnnData)
    source_per_mod_var_names = isinstance(source_var_names, dict)
    load_var_names = _get_var_names(
        adata, legacy_mudata_format=(not is_anndata and not source_per_mod_var_names)
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
            "`var_names` for the loaded `adata` does not match those of the `adata` used to train "
            "the model. For valid results, the former need to be the same and in the same order "
            "as the latter.",
            UserWarning,
            stacklevel=settings.warnings_stacklevel,
        )
