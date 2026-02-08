from __future__ import annotations

import logging
import warnings
from typing import TYPE_CHECKING
from uuid import uuid4

import h5py
import numpy as np
import pandas as pd
import scipy.sparse as sp_sparse
from anndata import AnnData
from anndata.abc import CSCDataset, CSRDataset
from anndata.io import read_elem
from mudata import MuData
from torch import as_tensor, sparse_csc_tensor, sparse_csr_tensor

from scvi import REGISTRY_KEYS, settings
from scvi.utils import is_package_installed

from . import _constants

if TYPE_CHECKING:
    from collections.abc import Iterator

    import numpy.typing as npt
    from pandas.api.types import CategoricalDtype
    from torch import Tensor

    from scvi._types import AnnOrMuData, MinifiedDataType

logger = logging.getLogger(__name__)

SparseDataset = (CSRDataset, CSCDataset)


def registry_key_to_default_dtype(key: str) -> type:
    """Returns the default dtype for a given registry key."""
    if key in [
        REGISTRY_KEYS.BATCH_KEY,
        REGISTRY_KEYS.LABELS_KEY,
        REGISTRY_KEYS.CAT_COVS_KEY,
        REGISTRY_KEYS.INDICES_KEY,
    ]:
        return np.int64

    return np.float32


def scipy_to_torch_sparse(x: sp_sparse.csr_matrix | sp_sparse.csc_matrix) -> Tensor:
    """Converts a SciPy sparse data structure to a sparse :class:`~torch.Tensor`.

    Parameters
    ----------
    x
        SciPy sparse data structure to convert. One of the following:

        * :class:`~scipy.sparse.csr_matrix`:
            Converted to a :class:`~torch.Tensor` constructed with
            :meth:`~torch.sparse_csr_tensor`.
        * :class:`~scipy.sparse.csc_matrix`:
            Converted to a :class:`~torch.Tensor` constructed with
            :meth:`~torch.sparse_csc_tensor`.

    Returns
    -------
    :class:`~torch.Tensor`
        A sparse tensor equivalent to `x`.
    """
    if isinstance(x, sp_sparse.csr_matrix):
        return sparse_csr_tensor(
            as_tensor(x.indptr),
            as_tensor(x.indices),
            as_tensor(x.data),
            size=x.shape,
        )
    elif isinstance(x, sp_sparse.csc_matrix):
        return sparse_csc_tensor(
            as_tensor(x.indptr),
            as_tensor(x.indices),
            as_tensor(x.data),
            size=x.shape,
        )
    else:
        raise TypeError(
            "`x` must be of type `scipy.sparse.csr_matrix` or `scipy.sparse.csc_matrix`."
        )


def get_anndata_attribute(
    adata: AnnOrMuData,
    attr_name: str,
    attr_key: str | None,
    mod_key: str | None = None,
) -> npt.NDArray | pd.DataFrame:
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
    elif isinstance(adata_attr, pd.DataFrame):
        if attr_key not in adata_attr.columns:
            raise ValueError(f"{attr_key} is not a valid column in adata.{attr_name}.")
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
    data: npt.NDArray | pd.DataFrame,
    attr_name: str,
    attr_key: str | None,
) -> None:
    """Sets the data in the AnnData object according to the attr_name and attr_key.

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


def _verify_and_correct_data_format(adata: AnnData, attr_name: str, attr_key: str | None):
    """Check data format and correct if necessary.

    Checks that the user's AnnData field is C_CONTIGUOUS and csr if it is dense numpy or sparse
    respectively.

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
        f"adata.{attr_name}[{attr_key}]" if attr_key is not None else f"adata.{attr_name}"
    )
    if sp_sparse.isspmatrix(data) and data.getformat() != "csr":
        warnings.warn(
            "Training will be faster when sparse matrix is formatted as CSR. It is safe to cast "
            "before model initialization.",
            UserWarning,
            stacklevel=settings.warnings_stacklevel,
        )
    elif isinstance(data, np.ndarray) and not data.flags["C_CONTIGUOUS"]:
        logger.debug(f"{data_loc_str} is not C_CONTIGUOUS. Overwriting to C_CONTIGUOUS.")
        data = np.asarray(data, order="C")
        _set_data_in_registry(adata, data, attr_name, attr_key)
    elif isinstance(data, pd.DataFrame) and not data.to_numpy().flags["C_CONTIGUOUS"]:
        logger.debug(f"{data_loc_str} is not C_CONTIGUOUS. Overwriting to C_CONTIGUOUS.")
        index = data.index
        vals = data.to_numpy()
        columns = data.columns
        data = pd.DataFrame(np.ascontiguousarray(vals), index=index, columns=columns)
        _set_data_in_registry(adata, data, attr_name, attr_key)


def _make_column_categorical(
    df: pd.DataFrame,
    column_key: str,
    alternate_column_key: str,
    categorical_dtype: str | CategoricalDtype | None = None,
):
    """Makes the data in column_key in DataFrame all categorical.

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
            f'Making .obs["{column_key}"] categorical failed. Expected categories: {mapping}. '
            f"Received categories: {received_categories}. "
        )
    df[alternate_column_key] = codes

    # make sure each category contains enough cells
    if np.min(counts) < 3:
        category = unique[np.argmin(counts)]
        warnings.warn(
            f"Category {category} in adata.obs['{alternate_column_key}'] has fewer than 3 cells. "
            "Models may not train properly.",
            UserWarning,
            stacklevel=settings.warnings_stacklevel,
        )

    return mapping


def _assign_adata_uuid(adata: AnnOrMuData, overwrite: bool = False) -> None:
    """Assigns a UUID unique to the AnnData object.

    If already present, the UUID is left alone, unless ``overwrite == True``.
    """
    if _constants._SCVI_UUID_KEY not in adata.uns or overwrite:
        adata.uns[_constants._SCVI_UUID_KEY] = str(uuid4())


def _check_nonnegative_integers(
    data: pd.DataFrame | npt.NDArray | sp_sparse.spmatrix | h5py.Dataset,
    n_to_check: int = 20,
):
    """Approximately checks values of data to ensure it is count data."""
    # for backed anndata
    if isinstance(data, h5py.Dataset) or isinstance(data, SparseDataset):
        data = data[:100]
    elif is_package_installed("dask"):
        import dask.array as da

        if isinstance(data, da.Array):
            data = data[:100, :100].compute()

    if isinstance(data, np.ndarray):
        data = data
    elif issubclass(type(data), sp_sparse.spmatrix):
        data = data.data
    elif isinstance(data, pd.DataFrame):
        data = data.to_numpy()
    else:
        raise TypeError("data type not understood")

    ret = True
    if data.shape[0] != 0:
        inds = np.random.choice(data.shape[0], size=(n_to_check,))
        check = data[inds]
        negative = np.any(check < 0)
        non_integer = np.any(check % 1 != 0)
        ret = not (negative or non_integer)
    return ret


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
        raise AssertionError("Cannot call ``_check_mudata_fully_paired`` with AnnData object.")
    for mod_key in mdata.mod:
        if not mdata.obsm[mod_key].all():
            raise ValueError(
                f"Detected unpaired observations in modality {mod_key}. "
                "Please make sure that data is fully paired in all MuData inputs. "
                "Either pad the unpaired modalities or take the intersection with "
                "muon.pp.intersect_obs()."
            )


def _get_adata_minify_type(adata: AnnData) -> MinifiedDataType | None:
    return adata.uns.get(_constants._ADATA_MINIFY_TYPE_UNS_KEY, None)


def _is_minified(adata: AnnOrMuData | str) -> bool:
    uns_key = _constants._ADATA_MINIFY_TYPE_UNS_KEY
    if isinstance(adata, AnnData):
        return adata.uns.get(uns_key, None) is not None
    elif isinstance(adata, MuData):
        return adata.uns.get(uns_key, None) is not None
    elif isinstance(adata, str):
        with h5py.File(adata) as fp:
            return uns_key in read_elem(fp["uns"]).keys()
    else:
        raise TypeError(f"Unsupported type: {type(adata)}")


def _check_fragment_counts(
    data: pd.DataFrame | npt.NDArray | sp_sparse.spmatrix | h5py.Dataset,
    n_to_check: int = 100,
):
    """Approximately checks values of data to ensure it is fragment count data."""
    # for backed anndata
    if isinstance(data, h5py.Dataset) or isinstance(data, SparseDataset):
        if data.shape[0] >= 400:
            data = data[:400]
        else:
            data = data[:]
    elif is_package_installed("dask"):
        import dask.array as da

        if isinstance(data, da.Array):
            if data.shape[0] >= 400:
                data = data[:400].compute()
            else:
                data = data[:].compute()

    # check that n_obs is greater than n_to_check
    if data.shape[0] < n_to_check:
        raise ValueError(
            f"adata.obs must have at least {n_to_check} observations. Consider reducing "
            "n_to_check."
        )

    inds = np.random.choice(data.shape[0], size=(n_to_check,))
    if isinstance(data, np.ndarray):
        data = data[inds]
    elif sp_sparse.issparse(data):
        data = data[inds].data
    elif isinstance(data, pd.DataFrame):
        data = data.iloc[inds].to_numpy()
    else:
        raise TypeError("data type not understood")

    binary = np.logical_not(np.any(data > 1))  # True if all values are 0 or 1
    non_fragments = np.count_nonzero(data == 1) < np.count_nonzero(
        data == 2
    )  # True if there are more 2s than 1s
    ret = not (non_fragments or binary)
    return ret


def _validate_adata_dataloader_input(
    model,
    adata: AnnOrMuData | None = None,
    dataloader: Iterator[dict[str, Tensor | None]] | None = None,
):
    """Validate that model uses adata or custom dataloader"""
    if adata is not None and dataloader is not None:
        raise ValueError("Only one of `adata` or `dataloader` can be provided.")
    elif (
        hasattr(model, "registry")
        and "setup_method_name" in model.registry.keys()
        and model.registry["setup_method_name"] == "setup_datamodule"
        and dataloader is None
    ):
        raise ValueError("`dataloader` must be provided.")
    return
