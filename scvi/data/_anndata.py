import logging
import os
import pickle
import sys
import warnings
from typing import Dict, List, Optional, Tuple, Union

import anndata
import numpy as np
import pandas as pd
import rich
from pandas.api.types import CategoricalDtype
from rich.console import Console
from scipy.sparse import isspmatrix

import scvi
from scvi import _CONSTANTS
from scvi._compat import Literal

from ._utils import (
    _check_nonnegative_integers,
    _compute_library_size_batch,
    _get_batch_mask_protein_data,
)

logger = logging.getLogger(__name__)


def get_from_registry(adata: anndata.AnnData, key: str) -> np.ndarray:
    """
    Returns the object in AnnData associated with the key in ``.uns['_scvi']['data_registry']``.

    Parameters
    ----------
    adata
        anndata object already setup with `scvi.data.setup_anndata()`
    key
        key of object to get from ``adata.uns['_scvi]['data_registry']``

    Returns
    -------
    The requested data

    Examples
    --------
    >>> import scvi
    >>> adata = scvi.data.cortex()
    >>> adata.uns['_scvi']['data_registry']
    {'X': ['_X', None],
    'batch_indices': ['obs', 'batch'],
    'local_l_mean': ['obs', '_scvi_local_l_mean'],
    'local_l_var': ['obs', '_scvi_local_l_var'],
    'labels': ['obs', 'labels']}
    >>> batch = get_from_registry(adata, "batch_indices")
    >>> batch
    array([[0],
           [0],
           [0],
           ...,
           [0],
           [0],
           [0]])
    """
    data_loc = adata.uns["_scvi"]["data_registry"][key]
    attr_name, attr_key = data_loc["attr_name"], data_loc["attr_key"]

    data = getattr(adata, attr_name)
    if attr_key != "None":
        if isinstance(data, pd.DataFrame):
            data = data.loc[:, attr_key]
        else:
            data = data[attr_key]
    if isinstance(data, pd.Series):
        data = data.to_numpy().reshape(-1, 1)
    return data


def setup_anndata(
    adata: anndata.AnnData,
    batch_key: Optional[str] = None,
    labels_key: Optional[str] = None,
    layer: Optional[str] = None,
    protein_expression_obsm_key: Optional[str] = None,
    protein_names_uns_key: Optional[str] = None,
    categorical_covariate_keys: Optional[List[str]] = None,
    continuous_covariate_keys: Optional[List[str]] = None,
    copy: bool = False,
) -> Optional[anndata.AnnData]:
    """
    Sets up :class:`~anndata.AnnData` object for `scvi` models.

    A mapping will be created between data fields used by `scvi` to their respective locations in adata.
    This method will also compute the log mean and log variance per batch for the library size prior.

    None of the data in adata are modified. Only adds fields to adata.

    Parameters
    ----------
    adata
        AnnData object containing raw counts. Rows represent cells, columns represent features.
    batch_key
        key in `adata.obs` for batch information. Categories will automatically be converted into integer
        categories and saved to `adata.obs['_scvi_batch']`. If `None`, assigns the same batch to all the data.
    labels_key
        key in `adata.obs` for label information. Categories will automatically be converted into integer
        categories and saved to `adata.obs['_scvi_labels']`. If `None`, assigns the same label to all the data.
    layer
        if not `None`, uses this as the key in `adata.layers` for raw count data.
    protein_expression_obsm_key
        key in `adata.obsm` for protein expression data, Required for :class:`~scvi.model.TOTALVI`.
    protein_names_uns_key
        key in `adata.uns` for protein names. If None, will use the column names of `adata.obsm[protein_expression_obsm_key]`
        if it is a DataFrame, else will assign sequential names to proteins. Only relevant but not required for :class:`~scvi.model.TOTALVI`.
    categorical_covariate_keys
        keys in `adata.obs` that correspond to categorical data. Used in some `scvi` models.
    continuous_covariate_keys
        keys in `adata.obs` that correspond to continuous data. Used in some `scvi` models.
    copy
        if `True`, a copy of adata is returned.

    Returns
    -------
    If ``copy``,  will return :class:`~anndata.AnnData`.
    Adds the following fields to adata:

    .uns['_scvi']
        `scvi` setup dictionary
    .obs['_local_l_mean']
        per batch library size mean
    .obs['_local_l_var']
        per batch library size variance
    .obs['_scvi_labels']
        labels encoded as integers
    .obs['_scvi_batch']
        batch encoded as integers

    Examples
    --------
    Example setting up a scanpy dataset with random gene data and no batch nor label information

    >>> import scanpy as sc
    >>> import scvi
    >>> import numpy as np
    >>> adata = scvi.data.synthetic_iid(run_setup_anndata=False)
    >>> adata
    AnnData object with n_obs × n_vars = 400 × 100
        obs: 'batch', 'labels'
        uns: 'protein_names'
        obsm: 'protein_expression'

    Filter cells and run preprocessing before `setup_anndata`

    >>> sc.pp.filter_cells(adata, min_counts = 0)

    Since no batch_key nor labels_key was passed, setup_anndata() will assume all cells have the same batch and label

    >>> scvi.data.setup_anndata(adata)
    INFO      No batch_key inputted, assuming all cells are same batch
    INFO      No label_key inputted, assuming all cells have same label
    INFO      Using data from adata.X
    INFO      Computing library size prior per batch
    INFO      Registered keys:['X', 'batch_indices', 'local_l_mean', 'local_l_var', 'labels']
    INFO      Successfully registered anndata object containing 400 cells, 100 vars, 1 batches, 1 labels, and 0 proteins. Also registered 0 extra categorical covariates and 0 extra continuous covariates.

    Example setting up scanpy dataset with random gene data, batch, and protein expression

    >>> adata = scvi.data.synthetic_iid(run_setup_anndata=False)
    >>> scvi.data.setup_anndata(adata, batch_key='batch', protein_expression_obsm_key='protein_expression')
    INFO      Using batches from adata.obs["batch"]
    INFO      No label_key inputted, assuming all cells have same label
    INFO      Using data from adata.X
    INFO      Computing library size prior per batch
    INFO      Using protein expression from adata.obsm['protein_expression']
    INFO      Generating sequential protein names
    INFO      Registered keys:['X', 'batch_indices', 'local_l_mean', 'local_l_var', 'labels', 'protein_expression']
    INFO      Successfully registered anndata object containing 400 cells, 100 vars, 2 batches, 1 labels, and 100 proteins. Also registered 0 extra categorical covariates and 0 extra continuous covariates.
    """
    if copy:
        adata = adata.copy()

    if adata.is_view:
        raise ValueError(
            "Please run `adata = adata.copy()` or use the copy option in this function."
        )

    adata.uns["_scvi"] = {}
    adata.uns["_scvi"]["scvi_version"] = scvi.__version__

    batch_key = _setup_batch(adata, batch_key)
    labels_key = _setup_labels(adata, labels_key)
    x_loc, x_key = _setup_x(adata, layer)
    local_l_mean_key, local_l_var_key = _setup_library_size(adata, batch_key, layer)

    data_registry = {
        _CONSTANTS.X_KEY: {"attr_name": x_loc, "attr_key": x_key},
        _CONSTANTS.BATCH_KEY: {"attr_name": "obs", "attr_key": batch_key},
        _CONSTANTS.LOCAL_L_MEAN_KEY: {"attr_name": "obs", "attr_key": local_l_mean_key},
        _CONSTANTS.LOCAL_L_VAR_KEY: {"attr_name": "obs", "attr_key": local_l_var_key},
        _CONSTANTS.LABELS_KEY: {"attr_name": "obs", "attr_key": labels_key},
    }

    if protein_expression_obsm_key is not None:
        protein_expression_obsm_key = _setup_protein_expression(
            adata, protein_expression_obsm_key, protein_names_uns_key, batch_key
        )
        data_registry[_CONSTANTS.PROTEIN_EXP_KEY] = {
            "attr_name": "obsm",
            "attr_key": protein_expression_obsm_key,
        }

    if categorical_covariate_keys is not None:
        cat_loc, cat_key = _setup_extra_categorical_covs(
            adata, categorical_covariate_keys
        )
        data_registry[_CONSTANTS.CAT_COVS_KEY] = {
            "attr_name": cat_loc,
            "attr_key": cat_key,
        }

    if continuous_covariate_keys is not None:
        cont_loc, cont_key = _setup_extra_continuous_covs(
            adata, continuous_covariate_keys
        )
        data_registry[_CONSTANTS.CONT_COVS_KEY] = {
            "attr_name": cont_loc,
            "attr_key": cont_key,
        }

    # add the data_registry to anndata
    _register_anndata(adata, data_registry_dict=data_registry)
    logger.debug("Registered keys:{}".format(list(data_registry.keys())))
    _setup_summary_stats(
        adata,
        batch_key,
        labels_key,
        protein_expression_obsm_key,
        categorical_covariate_keys,
        continuous_covariate_keys,
    )

    logger.info("Please do not further modify adata until model is trained.")

    _verify_and_correct_data_format(adata, data_registry)

    if copy:
        return adata


def _set_data_in_registry(adata, data, key):
    """
    Sets the data associated with key in adata.uns['_scvi']['data_registry'].keys() to data.

    Note: This is a dangerous method and will change the underlying data of the user's anndata
    Currently used to make the user's anndata C_CONTIGUOUS and csr if it is dense numpy
    or sparse respectively.

    Parameters
    ----------
    adata
        anndata object to change data of
    data
        data to change to
    key
        key in adata.uns['_scvi]['data_registry'].keys() associated with the data
    """
    data_loc = adata.uns["_scvi"]["data_registry"][key]
    attr_name, attr_key = data_loc["attr_name"], data_loc["attr_key"]

    if attr_key == "None":
        setattr(adata, attr_name, data)

    elif attr_key != "None":
        attribute = getattr(adata, attr_name)
        if isinstance(attribute, pd.DataFrame):
            attribute.loc[:, attr_key] = data
        else:
            attribute[attr_key] = data
        setattr(adata, attr_name, attribute)


def _verify_and_correct_data_format(adata, data_registry):
    """
    Will make sure that the user's anndata is C_CONTIGUOUS and csr if it is dense numpy or sparse respectively.

    Will iterate through all the keys of data_registry.

    Parameters
    ----------
    adata
        anndata to check
    data_registry
        data registry of anndata
    """
    keys_to_check = [_CONSTANTS.X_KEY, _CONSTANTS.PROTEIN_EXP_KEY]
    keys = [key for key in keys_to_check if key in data_registry.keys()]

    for k in keys:
        data = get_from_registry(adata, k)
        if isspmatrix(data) and (data.getformat() != "csr"):
            logger.warning(
                "Training will be faster when sparse matrix is formatted as CSR. It is safe to cast before model initialization."
            )
        elif isinstance(data, np.ndarray) and (data.flags["C_CONTIGUOUS"] is False):
            logger.debug(
                "{} is not C_CONTIGUOUS. Overwriting to C_CONTIGUOUS.".format(k)
            )
            data = np.asarray(data, order="C")
            _set_data_in_registry(adata, data, k)
        elif isinstance(data, pd.DataFrame) and (
            data.to_numpy().flags["C_CONTIGUOUS"] is False
        ):
            logger.debug(
                "{} is not C_CONTIGUOUS. Overwriting to C_CONTIGUOUS.".format(k)
            )
            index = data.index
            vals = data.to_numpy()
            columns = data.columns
            data = pd.DataFrame(
                np.ascontiguousarray(vals), index=index, columns=columns
            )
            _set_data_in_registry(adata, data, k)


def register_tensor_from_anndata(
    adata: anndata.AnnData,
    registry_key: str,
    adata_attr_name: Literal["obs", "var", "obsm", "varm", "uns"],
    adata_key_name: str,
    is_categorical: Optional[bool] = False,
    adata_alternate_key_name: Optional[str] = None,
):
    """
    Add another tensor to scvi data registry.

    This function is intended for contributors testing out new models.

    Parameters
    ----------
    adata
        AnnData with "_scvi" key in `.uns`
    registry_key
        Key for tensor in registry, which will be the key in the dataloader output
    adata_attr_name
        AnnData attribute with tensor
    adata_key_name
        key in adata_attr_name with data
    is_categorical
        Whether or not data is categorical
    adata_alternate_key_name
        Added key in adata_attr_name for categorical codes if `is_categorical` is True
    """
    if is_categorical is True:
        if adata_attr_name != "obs":
            raise ValueError("categorical handling only implemented for data in `.obs`")

    if is_categorical is True and adata_attr_name == "obs":
        adata_key_name = _make_obs_column_categorical(
            adata,
            column_key=adata_key_name,
            alternate_column_key=adata_alternate_key_name,
        )
    new_dict = {
        registry_key: {"attr_name": adata_attr_name, "attr_key": adata_key_name}
    }

    data_registry = adata.uns["_scvi"]["data_registry"]
    data_registry.update(new_dict)
    _verify_and_correct_data_format(adata, data_registry)


def transfer_anndata_setup(
    adata_source: Union[anndata.AnnData, dict],
    adata_target: anndata.AnnData,
    extend_categories: bool = False,
):
    """
    Transfer anndata setup from a source object to a target object.

    This handles encoding for categorical data and is useful in the case where an
    anndata object has been subsetted and a category is lost.

    Parameters
    ----------
    adata_source
        AnnData that has been setup with scvi. If `dict`, must be dictionary
        from source anndata containing scvi setup parameters.
    adata_target
        AnnData with equivalent organization as source, but possibly subsetted.
    extend_categories
        New categories in `adata_target` are added to the registry.
    """
    adata_target.uns["_scvi"] = {}

    if isinstance(adata_source, anndata.AnnData):
        _scvi_dict = adata_source.uns["_scvi"]
    else:
        _scvi_dict = adata_source
    data_registry = _scvi_dict["data_registry"]
    summary_stats = _scvi_dict["summary_stats"]

    # transfer version
    adata_target.uns["_scvi"]["scvi_version"] = _scvi_dict["scvi_version"]

    x_loc = data_registry[_CONSTANTS.X_KEY]["attr_name"]
    if x_loc == "layers":
        layer = data_registry[_CONSTANTS.X_KEY]["attr_key"]
    else:
        layer = None

    target_n_vars = adata_target.shape[1]

    if target_n_vars != summary_stats["n_vars"]:
        raise ValueError(
            "Number of vars in adata_target not the same as source. "
            + "Expected: {} Received: {}".format(target_n_vars, summary_stats["n_vars"])
        )

    # transfer batch and labels
    categorical_mappings = _scvi_dict["categorical_mappings"]
    _transfer_batch_and_labels(adata_target, categorical_mappings, extend_categories)

    batch_key = "_scvi_batch"
    labels_key = "_scvi_labels"

    # transfer protein_expression
    protein_expression_obsm_key = _transfer_protein_expression(
        _scvi_dict, adata_target, batch_key
    )

    # transfer X
    x_loc, x_key = _setup_x(adata_target, layer)
    local_l_mean_key, local_l_var_key = _setup_library_size(
        adata_target, batch_key, layer
    )
    target_data_registry = data_registry.copy()
    target_data_registry.update(
        {_CONSTANTS.X_KEY: {"attr_name": x_loc, "attr_key": x_key}}
    )
    # transfer extra categorical covs
    has_cat_cov = True if _CONSTANTS.CAT_COVS_KEY in data_registry.keys() else False
    if has_cat_cov:
        source_cat_dict = _scvi_dict["extra_categoricals"]["mappings"].copy()
        # extend categories
        if extend_categories:
            for key, mapping in source_cat_dict:
                for c in np.unique(adata_target.obs[key]):
                    if c not in mapping:
                        mapping = np.concatenate([mapping, [c]])
                source_cat_dict[key] = mapping
        # use ["keys"] to maintain correct order
        cat_loc, cat_key = _setup_extra_categorical_covs(
            adata_target,
            _scvi_dict["extra_categoricals"]["keys"],
            category_dict=source_cat_dict,
        )
        target_data_registry.update(
            {_CONSTANTS.CAT_COVS_KEY: {"attr_name": cat_loc, "attr_key": cat_key}}
        )
    else:
        source_cat_dict = None

    # transfer extra continuous covs
    has_cont_cov = True if _CONSTANTS.CONT_COVS_KEY in data_registry.keys() else False
    if has_cont_cov:
        obs_keys_names = _scvi_dict["extra_continuous_keys"]
        cont_loc, cont_key = _setup_extra_continuous_covs(
            adata_target, list(obs_keys_names)
        )
        target_data_registry.update(
            {_CONSTANTS.CONT_COVS_KEY: {"attr_name": cont_loc, "attr_key": cont_key}}
        )
    else:
        obs_keys_names = None

    # add the data_registry to anndata
    _register_anndata(adata_target, data_registry_dict=target_data_registry)
    logger.info("Registered keys:{}".format(list(target_data_registry.keys())))
    _setup_summary_stats(
        adata_target,
        batch_key,
        labels_key,
        protein_expression_obsm_key,
        source_cat_dict,
        obs_keys_names,
    )
    _verify_and_correct_data_format(adata_target, data_registry)


def _transfer_batch_and_labels(adata_target, categorical_mappings, extend_categories):

    for key, val in categorical_mappings.items():
        original_key = val["original_key"]
        if (key == original_key) and (original_key not in adata_target.obs.keys()):
            # case where original key and key are equal
            # caused when no batch or label key were given
            # when anndata_source was setup
            logger.info(
                ".obs[{}] not found in target, assuming every cell is same category".format(
                    original_key
                )
            )
            adata_target.obs[original_key] = np.zeros(
                adata_target.shape[0], dtype=np.int64
            )
        elif (key != original_key) and (original_key not in adata_target.obs.keys()):
            raise KeyError(
                '.obs["{}"] was used to setup source, but not found in target.'.format(
                    original_key
                )
            )
        mapping = val["mapping"].copy()
        # extend mapping for new categories
        if extend_categories:
            for c in np.unique(adata_target.obs[original_key]):
                if c not in mapping:
                    mapping = np.concatenate([mapping, [c]])
        cat_dtype = CategoricalDtype(categories=mapping, ordered=True)
        _make_obs_column_categorical(
            adata_target, original_key, key, categorical_dtype=cat_dtype
        )


def _transfer_protein_expression(_scvi_dict, adata_target, batch_key):
    data_registry = _scvi_dict["data_registry"]
    summary_stats = _scvi_dict["summary_stats"]

    has_protein = True if _CONSTANTS.PROTEIN_EXP_KEY in data_registry.keys() else False
    if has_protein is True:
        prev_protein_obsm_key = data_registry[_CONSTANTS.PROTEIN_EXP_KEY]["attr_key"]
        if prev_protein_obsm_key not in adata_target.obsm.keys():
            raise KeyError(
                "Can't find {} in adata_target.obsm for protein expressions.".format(
                    prev_protein_obsm_key
                )
            )
        else:
            assert (
                summary_stats["n_proteins"]
                == adata_target.obsm[prev_protein_obsm_key].shape[1]
            )
            protein_expression_obsm_key = prev_protein_obsm_key

            adata_target.uns["_scvi"]["protein_names"] = _scvi_dict["protein_names"]
            # batch mask totalVI
            batch_mask = _get_batch_mask_protein_data(
                adata_target, protein_expression_obsm_key, batch_key
            )

            # check if it's actually needed
            if np.sum([~b[1] for b in batch_mask.items()]) > 0:
                logger.info("Found batches with missing protein expression")
                adata_target.uns["_scvi"]["totalvi_batch_mask"] = batch_mask
    else:
        protein_expression_obsm_key = None

    return protein_expression_obsm_key


def _assert_key_in_obs(adata, key):
    assert key in adata.obs.keys(), "{} is not a valid key for in adata.obs".format(key)


def _setup_labels(adata, labels_key):
    # checking labels
    if labels_key is None:
        logger.info("No label_key inputted, assuming all cells have same label")
        labels_key = "_scvi_labels"
        adata.obs[labels_key] = np.zeros(adata.shape[0], dtype=np.int64)
        alt_key = labels_key
    else:
        _assert_key_in_obs(adata, labels_key)
        logger.info('Using labels from adata.obs["{}"]'.format(labels_key))
        alt_key = "_scvi_labels"
    labels_key = _make_obs_column_categorical(
        adata, column_key=labels_key, alternate_column_key=alt_key
    )
    return labels_key


def _setup_batch(adata, batch_key):
    # checking batch
    if batch_key is None:
        logger.info("No batch_key inputted, assuming all cells are same batch")
        batch_key = "_scvi_batch"
        adata.obs[batch_key] = np.zeros(adata.shape[0], dtype=np.int64)
        alt_key = batch_key
    else:
        _assert_key_in_obs(adata, batch_key)
        logger.info('Using batches from adata.obs["{}"]'.format(batch_key))
        alt_key = "_scvi_batch"
    batch_key = _make_obs_column_categorical(
        adata, column_key=batch_key, alternate_column_key=alt_key
    )
    return batch_key


def _setup_extra_categorical_covs(
    adata: anndata.AnnData,
    categorical_covariate_keys: List[str],
    category_dict: Dict[str, List[str]] = None,
):
    """
    Setup obsm df for extra categorical covariates.

    Parameters
    ----------
    adata
        AnnData to setup
    categorical_covariate_keys
        List of keys in adata.obs with categorical data
    category_dict
        Optional dictionary with keys being keys of categorical data in obs
        and values being precomputed categories for each obs vector

    """
    for key in categorical_covariate_keys:
        _assert_key_in_obs(adata, key)

    cat_loc = "obsm"
    cat_key = "_scvi_extra_categoricals"

    adata.uns["_scvi"]["extra_categoricals"] = {}

    categories = {}
    df = pd.DataFrame(index=adata.obs_names)
    for key in categorical_covariate_keys:
        if category_dict is None:
            categorical_obs = adata.obs[key].astype("category")
            mapping = categorical_obs.cat.categories.to_numpy(copy=True)
            categories[key] = mapping
        else:
            possible_cats = category_dict[key]
            categorical_obs = adata.obs[key].astype(
                CategoricalDtype(categories=possible_cats)
            )
        codes = categorical_obs.cat.codes
        df[key] = codes

    adata.obsm[cat_key] = df

    store_cats = categories if category_dict is None else category_dict
    adata.uns["_scvi"]["extra_categoricals"]["mappings"] = store_cats
    # this preserves the order of the keys added to the df
    adata.uns["_scvi"]["extra_categoricals"]["keys"] = categorical_covariate_keys

    # how many cats per key, in the preserved order
    n_cats_per_key = []
    for k in categorical_covariate_keys:
        n_cats_per_key.append(len(store_cats[k]))
    adata.uns["_scvi"]["extra_categoricals"]["n_cats_per_key"] = n_cats_per_key

    return cat_loc, cat_key


def _setup_extra_continuous_covs(
    adata: anndata.AnnData, continuous_covariate_keys: List[str]
):
    """
    Setup obsm df for extra continuous covariates.

    Parameters
    ----------
    adata
        AnnData to setup
    continuous_covariate_keys
        List of keys in adata.obs with continuous data
    """
    for key in continuous_covariate_keys:
        _assert_key_in_obs(adata, key)

    cont_loc = "obsm"
    cont_key = "_scvi_extra_continuous"

    series = []
    for key in continuous_covariate_keys:
        s = adata.obs[key]
        series.append(s)

    adata.obsm[cont_key] = pd.concat(series, axis=1)
    adata.uns["_scvi"]["extra_continuous_keys"] = adata.obsm[
        cont_key
    ].columns.to_numpy()

    return cont_loc, cont_key


def _make_obs_column_categorical(
    adata, column_key, alternate_column_key, categorical_dtype=None
):
    """
    Makes the data in column_key in obs all categorical.

    If adata.obs[column_key] is not categorical, will categorize
    and save to .obs[alternate_column_key]
    """
    if categorical_dtype is None:
        categorical_obs = adata.obs[column_key].astype("category")
    else:
        categorical_obs = adata.obs[column_key].astype(categorical_dtype)

    # put codes in .obs[alternate_column_key]
    codes = categorical_obs.cat.codes
    mapping = categorical_obs.cat.categories.to_numpy(copy=True)
    if -1 in np.unique(codes):
        received_categories = adata.obs[column_key].astype("category").cat.categories
        raise ValueError(
            'Making .obs["{}"] categorical failed. Expected categories: {}. '
            "Received categories: {}. ".format(column_key, mapping, received_categories)
        )
    adata.obs[alternate_column_key] = codes

    # store categorical mappings
    store_dict = {
        alternate_column_key: {"original_key": column_key, "mapping": mapping}
    }
    if "categorical_mappings" not in adata.uns["_scvi"].keys():
        adata.uns["_scvi"].update({"categorical_mappings": store_dict})
    else:
        adata.uns["_scvi"]["categorical_mappings"].update(store_dict)

    # make sure each category contains enough cells
    unique, counts = np.unique(adata.obs[alternate_column_key], return_counts=True)
    if np.min(counts) < 3:
        category = unique[np.argmin(counts)]
        warnings.warn(
            "Category {} in adata.obs['{}'] has fewer than 3 cells. SCVI may not train properly.".format(
                category, alternate_column_key
            )
        )
    # possible check for continuous?
    if len(unique) > (adata.shape[0] / 3):
        warnings.warn(
            "Is adata.obs['{}'] continuous? SCVI doesn't support continuous obs yet."
        )
    return alternate_column_key


def _setup_protein_expression(
    adata, protein_expression_obsm_key, protein_names_uns_key, batch_key
):
    assert (
        protein_expression_obsm_key in adata.obsm.keys()
    ), "{} is not a valid key in adata.obsm".format(protein_expression_obsm_key)

    logger.info(
        "Using protein expression from adata.obsm['{}']".format(
            protein_expression_obsm_key
        )
    )
    pro_exp = adata.obsm[protein_expression_obsm_key]
    if _check_nonnegative_integers(pro_exp) is False:
        warnings.warn(
            "adata.obsm[{}] does not contain unnormalized count data. Are you sure this is what you want?".format(
                protein_expression_obsm_key
            )
        )
    # setup protein names
    if protein_names_uns_key is None and isinstance(
        adata.obsm[protein_expression_obsm_key], pd.DataFrame
    ):
        logger.info(
            "Using protein names from columns of adata.obsm['{}']".format(
                protein_expression_obsm_key
            )
        )
        protein_names = list(adata.obsm[protein_expression_obsm_key].columns)
    elif protein_names_uns_key is not None:
        logger.info(
            "Using protein names from adata.uns['{}']".format(protein_names_uns_key)
        )
        protein_names = adata.uns[protein_names_uns_key]
    else:
        logger.info("Generating sequential protein names")
        protein_names = np.arange(adata.obsm[protein_expression_obsm_key].shape[1])

    adata.uns["_scvi"]["protein_names"] = protein_names

    # batch mask totalVI
    batch_mask = _get_batch_mask_protein_data(
        adata, protein_expression_obsm_key, batch_key
    )

    # check if it's actually needed
    if np.sum([~b[1] for b in batch_mask.items()]) > 0:
        logger.info("Found batches with missing protein expression")
        adata.uns["_scvi"]["totalvi_batch_mask"] = batch_mask
    return protein_expression_obsm_key


def _setup_x(adata, layer):
    if layer is not None:
        assert (
            layer in adata.layers.keys()
        ), "{} is not a valid key in adata.layers".format(layer)
        logger.info('Using data from adata.layers["{}"]'.format(layer))
        x_loc = "layers"
        x_key = layer
        x = adata.layers[x_key]
    else:
        logger.info("Using data from adata.X")
        x_loc = "X"
        x_key = "None"
        x = adata.X

    if _check_nonnegative_integers(x) is False:
        logger_data_loc = (
            "adata.X" if layer is None else "adata.layers[{}]".format(layer)
        )
        warnings.warn(
            "{} does not contain unnormalized count data. Are you sure this is what you want?".format(
                logger_data_loc
            )
        )

    return x_loc, x_key


def _setup_library_size(adata, batch_key, layer):
    # computes the library size per batch
    logger.info("Computing library size prior per batch")
    local_l_mean_key = "_scvi_local_l_mean"
    local_l_var_key = "_scvi_local_l_var"
    _compute_library_size_batch(
        adata,
        batch_key=batch_key,
        local_l_mean_key=local_l_mean_key,
        local_l_var_key=local_l_var_key,
        layer=layer,
    )
    return local_l_mean_key, local_l_var_key


def _setup_summary_stats(
    adata,
    batch_key,
    labels_key,
    protein_expression_obsm_key,
    categorical_covariate_keys,
    continuous_covariate_keys,
):
    categorical_mappings = adata.uns["_scvi"]["categorical_mappings"]
    n_batch = len(np.unique(categorical_mappings[batch_key]["mapping"]))
    n_cells = adata.shape[0]
    n_vars = adata.shape[1]
    n_labels = len(np.unique(categorical_mappings[labels_key]["mapping"]))

    if protein_expression_obsm_key is not None:
        n_proteins = adata.obsm[protein_expression_obsm_key].shape[1]
    else:
        n_proteins = 0

    if categorical_covariate_keys is not None:
        n_cat_covs = len(categorical_covariate_keys)
    else:
        n_cat_covs = 0

    if continuous_covariate_keys is not None:
        n_cont_covs = len(continuous_covariate_keys)
    else:
        n_cont_covs = 0

    summary_stats = {
        "n_batch": n_batch,
        "n_cells": n_cells,
        "n_vars": n_vars,
        "n_labels": n_labels,
        "n_proteins": n_proteins,
        "n_continuous_covs": n_cont_covs,
    }
    adata.uns["_scvi"]["summary_stats"] = summary_stats
    logger.info(
        "Successfully registered anndata object containing {} cells, {} vars, "
        "{} batches, {} labels, and {} proteins. Also registered {} extra categorical "
        "covariates and {} extra continuous covariates.".format(
            n_cells, n_vars, n_batch, n_labels, n_proteins, n_cat_covs, n_cont_covs
        )
    )

    return summary_stats


def _register_anndata(adata, data_registry_dict: Dict[str, Tuple[str, str]]):
    """
    Registers the AnnData object by adding data_registry_dict to adata.uns['_scvi']['data_registry'].

    Format of data_registry_dict is: {<scvi_key>: (<anndata dataframe>, <dataframe key> )}

    Parameters
    ----------
    adata
        anndata object
    data_registry_dict
        dictionary mapping keys used by scvi.model to their respective location in adata.

    Examples
    --------
    >>> data_dict = {"batch" :("obs", "batch_idx"), "X": ("_X", None)}
    >>> _register_anndata(adata, data_dict)
    """
    adata.uns["_scvi"]["data_registry"] = data_registry_dict.copy()


def view_anndata_setup(source: Union[anndata.AnnData, dict, str]):
    """
    Prints setup anndata.

    Parameters
    ----------
    source
        Either AnnData, path to saved AnnData, path to folder with adata.h5ad,
        or scvi-setup-dict (adata.uns['_scvi'])

    Examples
    --------
    >>> scvi.data.view_anndata_setup(adata)
    >>> scvi.data.view_anndata_setup('saved_model_folder/adata.h5ad')
    >>> scvi.data.view_anndata_setup('saved_model_folder/')
    >>> scvi.data.view_anndata_setup(adata.uns['_scvi'])
    """
    if isinstance(source, anndata.AnnData):
        adata = source
    elif isinstance(source, str):
        # check if user passed in folder or anndata
        if source.endswith("h5ad"):
            path = source
            adata = anndata.read(path)
        else:
            path = os.path.join(source, "adata.h5ad")
            if os.path.exists(path):
                adata = anndata.read(path)
            else:
                path = os.path.join(source, "attr.pkl")
                with open(path, "rb") as handle:
                    adata = None
                    setup_dict = pickle.load(handle)["scvi_setup_dict_"]
    elif isinstance(source, dict):
        adata = None
        setup_dict = source
    else:
        raise ValueError(
            "Invalid source passed in. Must be either AnnData, path to saved AnnData, "
            + "path to folder with adata.h5ad or scvi-setup-dict (adata.uns['_scvi'])"
        )

    if adata is not None:
        if "_scvi" not in adata.uns.keys():
            raise ValueError("Please run setup_anndata() on your adata first.")
        setup_dict = adata.uns["_scvi"]

    summary_stats = setup_dict["summary_stats"]
    data_registry = setup_dict["data_registry"]
    mappings = setup_dict["categorical_mappings"]
    version = setup_dict["scvi_version"]

    rich.print("Anndata setup with scvi-tools version {}.".format(version))

    n_cat = 0
    n_covs = 0
    if "extra_categorical_mappings" in setup_dict.keys():
        n_cat = len(setup_dict["extra_categoricals"]["mappings"])
    if "extra_continuous_keys" in setup_dict.keys():
        n_covs = len(setup_dict["extra_continuous_keys"])

    in_colab = "google.colab" in sys.modules
    force_jupyter = None if not in_colab else True
    console = Console(force_jupyter=force_jupyter)
    t = rich.table.Table(title="Data Summary")
    t.add_column(
        "Data", justify="center", style="dodger_blue1", no_wrap=True, overflow="fold"
    )
    t.add_column(
        "Count", justify="center", style="dark_violet", no_wrap=True, overflow="fold"
    )
    data_summary = {
        "Cells": summary_stats["n_cells"],
        "Vars": summary_stats["n_vars"],
        "Labels": summary_stats["n_labels"],
        "Batches": summary_stats["n_batch"],
        "Proteins": summary_stats["n_proteins"],
        "Extra Categorical Covariates": n_cat,
        "Extra Continuous Covariates": n_covs,
    }
    for data, count in data_summary.items():
        t.add_row(data, str(count))
    console.print(t)

    t = rich.table.Table(title="SCVI Data Registry")
    t.add_column(
        "Data", justify="center", style="dodger_blue1", no_wrap=True, overflow="fold"
    )
    t.add_column(
        "scvi-tools Location",
        justify="center",
        style="dark_violet",
        no_wrap=True,
        overflow="fold",
    )

    for scvi_data_key, data_loc in data_registry.items():
        attr_name = data_loc["attr_name"]
        attr_key = data_loc["attr_key"]
        if attr_key == "None":
            scvi_data_str = "adata.{}".format(attr_name)
        else:
            scvi_data_str = "adata.{}['{}']".format(attr_name, attr_key)

        t.add_row(scvi_data_key, scvi_data_str)

    console.print(t)

    t = _categorical_mappings_table("Label Categories", "_scvi_labels", mappings)
    console.print(t)
    t = _categorical_mappings_table("Batch Categories", "_scvi_batch", mappings)
    console.print(t)

    if "extra_categoricals" in setup_dict.keys():
        t = _extra_categoricals_table(setup_dict)
        console.print(t)

    if "extra_continuous_keys" in setup_dict.keys():
        t = _extra_continuous_table(adata, setup_dict)
        console.print(t)


def _extra_categoricals_table(setup_dict: dict):
    """Returns rich.table.Table with info on extra categorical variables."""
    t = rich.table.Table(title="Extra Categorical Variables")
    t.add_column(
        "Source Location",
        justify="center",
        style="dodger_blue1",
        no_wrap=True,
        overflow="fold",
    )
    t.add_column(
        "Categories", justify="center", style="green", no_wrap=True, overflow="fold"
    )
    t.add_column(
        "scvi-tools Encoding",
        justify="center",
        style="dark_violet",
        no_wrap=True,
        overflow="fold",
    )
    for key, mappings in setup_dict["extra_categoricals"]["mappings"].items():
        for i, mapping in enumerate(mappings):
            if i == 0:
                t.add_row("adata.obs['{}']".format(key), str(mapping), str(i))
            else:
                t.add_row("", str(mapping), str(i))
        t.add_row("", "")
    return t


def _extra_continuous_table(adata: Optional[anndata.AnnData], setup_dict: dict):
    """Returns rich.table.Table with info on extra continuous variables."""
    t = rich.table.Table(title="Extra Continuous Variables")
    t.add_column(
        "Source Location",
        justify="center",
        style="dodger_blue1",
        no_wrap=True,
        overflow="fold",
    )
    if adata is not None:
        t.add_column(
            "Range",
            justify="center",
            style="dark_violet",
            no_wrap=True,
            overflow="fold",
        )
        cont_covs = scvi.data.get_from_registry(adata, "cont_covs")
        for cov in cont_covs.iteritems():
            col_name, values = cov[0], cov[1]
            min_val = np.min(values)
            max_val = np.max(values)
            t.add_row(
                "adata.obs['{}']".format(col_name),
                "{:.20g} -> {:.20g}".format(min_val, max_val),
            )
    else:
        for key in setup_dict["extra_continuous_keys"]:
            t.add_row("adata.obs['{}']".format(key))
    return t


def _categorical_mappings_table(title: str, scvi_column: str, mappings: dict):
    """
    Returns rich.table.Table with info on a categorical variable.

    Parameters
    ----------
    title
        title of table
    scvi_column
        column used by scvi for categorical representation
    mappings
        output of adata.uns['_scvi']['categorical_mappings'], containing mapping
        between scvi_column and original column and categories
    """
    source_key = mappings[scvi_column]["original_key"]
    mapping = mappings[scvi_column]["mapping"]
    t = rich.table.Table(title=title)
    t.add_column(
        "Source Location",
        justify="center",
        style="dodger_blue1",
        no_wrap=True,
        overflow="fold",
    )
    t.add_column(
        "Categories", justify="center", style="green", no_wrap=True, overflow="fold"
    )
    t.add_column(
        "scvi-tools Encoding",
        justify="center",
        style="dark_violet",
        no_wrap=True,
        overflow="fold",
    )
    for i, cat in enumerate(mapping):
        if i == 0:
            t.add_row("adata.obs['{}']".format(source_key), str(cat), str(i))
        else:
            t.add_row("", str(cat), str(i))
    return t


def _check_anndata_setup_equivalence(adata_source, adata_target):
    """
    Checks if target setup is equivalent to source.

    Parameters
    ----------
    adata_source
        Either AnnData already setup or scvi_setup_dict as the source
    adata_target
        Target AnnData to check setup equivalence
    """
    if isinstance(adata_source, anndata.AnnData):
        _scvi_dict = adata_source.uns["_scvi"]
    else:
        _scvi_dict = adata_source
    adata = adata_target

    stats = _scvi_dict["summary_stats"]

    target_n_vars = adata.shape[1]
    error_msg = (
        "Number of {} in anndata different from initial anndata used for training."
    )
    if target_n_vars != stats["n_vars"]:
        raise ValueError(error_msg.format("vars"))

    error_msg = (
        "There are more {} categories in the data than were originally registered. "
        + "Please check your {} categories as well as adata.uns['_scvi']['categorical_mappings']."
    )
    self_categoricals = _scvi_dict["categorical_mappings"]
    self_batch_mapping = self_categoricals["_scvi_batch"]["mapping"]

    adata_categoricals = adata.uns["_scvi"]["categorical_mappings"]
    adata_batch_mapping = adata_categoricals["_scvi_batch"]["mapping"]

    # check if mappings are equal or needs transfer
    transfer_setup = _needs_transfer(self_batch_mapping, adata_batch_mapping, "batch")
    self_labels_mapping = self_categoricals["_scvi_labels"]["mapping"]
    adata_labels_mapping = adata_categoricals["_scvi_labels"]["mapping"]

    transfer_setup = transfer_setup or _needs_transfer(
        self_labels_mapping, adata_labels_mapping, "label"
    )

    # validate any extra categoricals
    if "extra_categoricals" in _scvi_dict.keys():
        target_dict = adata.uns["_scvi"]["extra_categoricals"]
        source_dict = _scvi_dict["extra_categoricals"]
        # check that order of keys setup is same
        if not np.array_equal(target_dict["keys"], source_dict["keys"]):
            error_msg = (
                "Registered categorical key order mismatch between "
                + "the anndata used to train and the anndata passed in."
                + "Expected categories & order {}. Received {}.\n"
            )
            raise ValueError(error_msg.format(source_dict["keys"], target_dict["keys"]))
        # check mappings are equivalent
        target_extra_cat_maps = adata.uns["_scvi"]["extra_categoricals"]["mappings"]
        for key, val in source_dict["mappings"].items():
            target_map = target_extra_cat_maps[key]
            transfer_setup = transfer_setup or _needs_transfer(val, target_map, key)
    # validate any extra continuous covs
    if "extra_continuous_keys" in _scvi_dict.keys():
        if "extra_continuous_keys" not in adata.uns["_scvi"].keys():
            raise ValueError('extra_continuous_keys not in adata.uns["_scvi"]')
        target_cont_keys = adata.uns["_scvi"]["extra_continuous_keys"]
        source_cont_keys = _scvi_dict["extra_continuous_keys"]
        n_keys = len(target_cont_keys)
        if np.sum(source_cont_keys == target_cont_keys) != n_keys:
            raise ValueError(
                "extra_continous_keys are not the same between source and target"
            )
    if transfer_setup:
        transfer_anndata_setup(adata_source, adata_target)


def _needs_transfer(mapping1, mapping2, category):
    needs_transfer = False
    error_msg = (
        "Categorial encoding for {} is not the same between "
        + "the anndata used to train and the anndata passed in. "
        + "Categorical encoding needs to be same elements, same order, and same datatype.\n"
        + "Expected categories: {}. Received categories: {}.\n"
    )
    warning_msg = (
        "Categorical encoding for {} is similar but not equal between "
        + "the anndata used to train and the anndata passed in. "
        + "Will attempt transfer. Expected categories: {}. Received categories: {}.\n "
    )
    if _is_equal_mapping(mapping1, mapping2):
        needs_transfer = False
    elif _is_similar_mapping(mapping1, mapping2):
        needs_transfer = True
        logger.warning(warning_msg.format(category, mapping1, mapping2))
    else:
        raise ValueError(error_msg.format(category, mapping1, mapping2))
    return needs_transfer


def _is_similar_mapping(mapping1, mapping2):
    """Returns True if mapping2 is a subset of mapping1."""
    if len(set(mapping2) - set(mapping1)) == 0:
        return True
    else:
        return False


def _is_equal_mapping(mapping1, mapping2):
    return pd.Index(mapping1).equals(pd.Index(mapping2))
