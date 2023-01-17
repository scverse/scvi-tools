import logging
import warnings
from typing import Optional, Union
from uuid import uuid4

import h5py
import jax
import jax.numpy as jnp
import numpy as np
import pandas as pd
import scipy.sparse as sp_sparse
from anndata import AnnData
from anndata._core.sparse_dataset import SparseDataset

# TODO use the experimental api once we lower bound to anndata 0.8
try:
    from anndata.experimental import read_elem
except ImportError:
    from anndata._io.specs import read_elem

from mudata import MuData
from pandas.api.types import CategoricalDtype

from scvi._types import AnnOrMuData, LatentDataType

from . import _constants

logger = logging.getLogger(__name__)


def get_anndata_attribute(
    adata: AnnOrMuData,
    attr_name: str,
    attr_key: Optional[str],
    mod_key: Optional[str] = None,
) -> Union[np.ndarray, pd.DataFrame]:
    """Returns the requested data from a given AnnData/MuData object."""
    if mod_key is not None:
        if isinstance(adata, AnnData):
            raise ValueError(f"Cannot access modality {mod_key} on an AnnData object.")
        if mod_key not in adata.mod:
            raise ValueError(f"{mod_key} is not a valid modality in adata.mod.")
        adata = adata.mod[mod_key]
    adata_attr = getattr(adata, attr_name)
    if attr_key is None:
        field = adata_attr
    else:
        if isinstance(adata_attr, pd.DataFrame):
            if attr_key not in adata_attr.columns:
                raise ValueError(
                    f"{attr_key} is not a valid column in adata.{attr_name}."
                )
            field = adata_attr.loc[:, attr_key]
        else:
            if attr_key not in adata_attr.keys():
                raise ValueError(f"{attr_key} is not a valid key in adata.{attr_name}.")
            field = adata_attr[attr_key]
    if isinstance(field, pd.Series):
        field = field.to_numpy().reshape(-1, 1)
    return field


def _set_data_in_registry(
    adata: AnnData,
    data: Union[np.ndarray, pd.DataFrame],
    attr_name: str,
    attr_key: Optional[str],
):
    """
    Sets the data in the AnnData object according to the attr_name and attr_key.

    Note: This is a dangerous method and will change the underlying data of the user's anndata
    Currently used to make the user's anndata C_CONTIGUOUS and csr if it is dense numpy
    or sparse respectively.

    Parameters
    ----------
    adata
        AnnData object to change data of.
    data
        Data to change to.
    attr_name
        Attribute name of AnnData object to store data in.
    attr_key
        Key in AnnData attribute under which to store data in.
    """
    if attr_key is None:
        setattr(adata, attr_name, data)

    elif attr_key is not None:
        attribute = getattr(adata, attr_name)
        if isinstance(attribute, pd.DataFrame):
            attribute.loc[:, attr_key] = data
        else:
            attribute[attr_key] = data
        setattr(adata, attr_name, attribute)


def _verify_and_correct_data_format(
    adata: AnnData, attr_name: str, attr_key: Optional[str]
):
    """
    Will make sure that the user's AnnData field is C_CONTIGUOUS and csr if it is dense numpy or sparse respectively.

    Parameters
    ----------
    adata
        AnnData object to check.
    attr_name
        Attribute name where data is stored.
    attr_key
        Attribute key where data is stored, if applicable.
    """
    data = get_anndata_attribute(adata, attr_name, attr_key)
    data_loc_str = (
        f"adata.{attr_name}[{attr_key}]"
        if attr_key is not None
        else f"adata.{attr_name}"
    )
    if sp_sparse.isspmatrix(data) and (data.getformat() != "csr"):
        warnings.warn(
            "Training will be faster when sparse matrix is formatted as CSR. It is safe to cast before model initialization."
        )
    elif isinstance(data, np.ndarray) and (data.flags["C_CONTIGUOUS"] is False):
        logger.debug(
            f"{data_loc_str} is not C_CONTIGUOUS. Overwriting to C_CONTIGUOUS."
        )
        data = np.asarray(data, order="C")
        _set_data_in_registry(adata, data, attr_name, attr_key)
    elif isinstance(data, pd.DataFrame) and (
        data.to_numpy().flags["C_CONTIGUOUS"] is False
    ):
        logger.debug(
            f"{data_loc_str} is not C_CONTIGUOUS. Overwriting to C_CONTIGUOUS."
        )
        index = data.index
        vals = data.to_numpy()
        columns = data.columns
        data = pd.DataFrame(np.ascontiguousarray(vals), index=index, columns=columns)
        _set_data_in_registry(adata, data, attr_name, attr_key)


def _make_column_categorical(
    df: pd.DataFrame,
    column_key: str,
    alternate_column_key: str,
    categorical_dtype: Optional[Union[str, CategoricalDtype]] = None,
):
    """
    Makes the data in column_key in DataFrame all categorical.

    Categorizes df[column_key], then saves category codes to
    df[alternate_column_key] and returns the category mappings.
    """
    if categorical_dtype is None:
        categorical_obs = df[column_key].astype("category")
    else:
        categorical_obs = df[column_key].astype(categorical_dtype)

    # put codes in .obs[alternate_column_key]
    codes = categorical_obs.cat.codes
    unique, counts = np.unique(codes, return_counts=True)
    mapping = categorical_obs.cat.categories.to_numpy(copy=True)
    if -1 in unique:
        received_categories = df[column_key].astype("category").cat.categories
        raise ValueError(
            'Making .obs["{}"] categorical failed. Expected categories: {}. '
            "Received categories: {}. ".format(column_key, mapping, received_categories)
        )
    df[alternate_column_key] = codes

    # make sure each category contains enough cells
    if np.min(counts) < 3:
        category = unique[np.argmin(counts)]
        warnings.warn(
            "Category {} in adata.obs['{}'] has fewer than 3 cells. Models may not train properly.".format(
                category, alternate_column_key
            )
        )

    return mapping


def _assign_adata_uuid(adata: AnnOrMuData, overwrite: bool = False) -> None:
    """
    Assigns a UUID unique to the AnnData object.

    If already present, the UUID is left alone, unless ``overwrite == True``.
    """
    if _constants._SCVI_UUID_KEY not in adata.uns or overwrite:
        adata.uns[_constants._SCVI_UUID_KEY] = str(uuid4())


def _check_nonnegative_integers(
    data: Union[pd.DataFrame, np.ndarray, sp_sparse.spmatrix, h5py.Dataset],
    n_to_check: int = 20,
):
    """Approximately checks values of data to ensure it is count data."""
    # for backed anndata
    if isinstance(data, h5py.Dataset) or isinstance(data, SparseDataset):
        data = data[:100]

    if isinstance(data, np.ndarray):
        data = data
    elif issubclass(type(data), sp_sparse.spmatrix):
        data = data.data
    elif isinstance(data, pd.DataFrame):
        data = data.to_numpy()
    else:
        raise TypeError("data type not understood")

    ret = True
    if len(data) != 0:
        inds = np.random.choice(len(data), size=(n_to_check,))
        check = jax.device_put(data.flat[inds], device=jax.devices("cpu")[0])
        negative, non_integer = _is_not_count_val(check)
        ret = not (negative or non_integer)
    return ret


@jax.jit
def _is_not_count_val(data: jnp.ndarray):
    negative = jnp.any(data < 0)
    non_integer = jnp.any(data % 1 != 0)

    return negative, non_integer


def _check_if_view(adata: AnnOrMuData, copy_if_view: bool = False):
    if adata.is_view:
        if copy_if_view:
            logger.info("Received view of anndata, making copy.")
            adata._init_as_actual(adata.copy())
            # Reassign AnnData UUID to produce a separate AnnDataManager.
            _assign_adata_uuid(adata, overwrite=True)
        else:
            raise ValueError("Please run `adata = adata.copy()`")
    elif isinstance(adata, MuData):
        for mod_key in adata.mod.keys():
            mod_adata = adata.mod[mod_key]
            _check_if_view(mod_adata, copy_if_view)


def _check_mudata_fully_paired(mdata: MuData):
    if isinstance(mdata, AnnData):
        raise AssertionError(
            "Cannot call ``_check_mudata_fully_paired`` with AnnData object."
        )
    for mod_key in mdata.mod:
        if not mdata.obsm[mod_key].all():
            raise ValueError(
                f"Detected unpaired observations in modality {mod_key}. "
                "Please make sure that data is fully paired in all MuData inputs. "
                "Either pad the unpaired modalities or take the intersection with muon.pp.intersect_obs()."
            )


def _get_latent_adata_type(adata: AnnData) -> Union[LatentDataType, None]:
    return adata.uns.get(_constants._ADATA_LATENT_UNS_KEY, None)


def _is_latent(adata: Union[AnnData, str]) -> bool:
    uns_key = _constants._ADATA_LATENT_UNS_KEY
    if isinstance(adata, AnnData):
        return adata.uns.get(uns_key, None) is not None
    elif isinstance(adata, str):
        with h5py.File(adata) as fp:
            return uns_key in read_elem(fp["uns"]).keys()
    else:
        raise TypeError(f"Unsupported type: {type(adata)}")
