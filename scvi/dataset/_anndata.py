import copy
import warnings
import numpy as np
import logging
import pandas as pd
import anndata

from typing import Dict, Tuple, Optional
from scvi.dataset._anndata_utils import (
    _compute_library_size_batch,
    _check_nonnegative_integers,
    _get_batch_mask_protein_data,
)
from scvi import _CONSTANTS_

logger = logging.getLogger(__name__)


def _register_anndata(adata, data_registry_dict: Dict[str, Tuple[str, str]]):
    """Registers the AnnData object by adding data_registry_dict to adata.uns

    Format is: {<scvi_key>: (<anndata dataframe>, <dataframe key> )}
    Example:
    {"batch" :("obs", "batch_idx")}
    {"X": ("_X", None)}

    Parameters
    ----------
    adata
        anndata object
    data_registry_dict
        dictionary mapping keys used by scvi models to their respective location in adata.

    """
    for df, df_key in data_registry_dict.values():
        if df_key is not None:
            assert df_key in getattr(
                adata, df
            ), "anndata.{} has no attribute '{}'".format(df, df_key)
        else:
            assert hasattr(adata, df) is True, "anndata has no attribute '{}'".format(
                df_key
            )
    adata.uns["scvi_data_registry"] = copy.copy(data_registry_dict)


def get_from_registry(adata, key: str):
    """Returns an the object in Anndata associated the key in adata.uns['scvi_data_registry']

    Parameters
    ----------
    adata
        anndata object
    key
        key of object to get from adata.uns['scvi_data_registry']

    Returns
    -------

    """
    data_loc = adata.uns["scvi_data_registry"][key]
    df, df_key = data_loc[0], data_loc[1]

    if df_key == "":
        df_key = None
    data = getattr(adata, df)
    if df_key is not None:
        data = data[df_key]
    if isinstance(data, pd.Series):
        data = np.array(data.values).reshape(adata.shape[0], -1)
    return data


def setup_anndata(
    adata,
    batch_key: str = None,
    labels_key: str = None,
    X_layers_key: str = None,
    protein_expression_obsm_key: str = None,
    protein_names_uns_key: str = None,
    copy: bool = False,
) -> Optional[anndata.AnnData]:
    """Sets up anndata object for scVI models.

    This method will compute the log mean and log variance per batch.
    A mapping will be created between in

    Parameters
    ----------
    adata
        anndata object containing raw counts
    batch_key
        key in adata.obs for batch information. Will automatically be converted into integer categories
    labels_key
        key in adata.obs for label information. Will automatically be converted into integer categories
    X_layers_key
        if not None, uses this as the key in adata.layers for raw count
    protein_expression_obsm_key
        key in adata.obsm for protein expression data
    protein_names_uns_key
        key in adata.uns for protein names
    copy
        if True, a copy of anndata is returned

    Returns
    -------
    """
    if adata.is_view:
        raise ValueError("adata cannot be a view of an AnnData object.")

    if copy:
        adata = adata.copy()

    batch_key = _setup_batch(adata, batch_key)
    labels_key = _setup_labels(adata, labels_key)
    X_loc, X_key = _setup_X(adata, X_layers_key)
    local_l_mean_key, local_l_var_key = _setup_library_size(
        adata, batch_key, X_layers_key
    )

    data_registry = {
        _CONSTANTS_.X_KEY: [X_loc, X_key],
        _CONSTANTS_.BATCH_KEY: ["obs", batch_key],
        _CONSTANTS_.LOCAL_L_MEAN_KEY: ["obs", local_l_mean_key],
        _CONSTANTS_.LOCAL_L_VAR_KEY: ["obs", local_l_var_key],
        _CONSTANTS_.LABELS_KEY: ["obs", labels_key],
    }

    if protein_expression_obsm_key is not None:
        protein_expression_obsm_key = _setup_protein_expression(
            adata, protein_expression_obsm_key, protein_names_uns_key, batch_key
        )
        data_registry[_CONSTANTS_.PROTEIN_EXP_KEY] = [
            "obsm",
            protein_expression_obsm_key,
        ]

    # add the data_registry to anndata
    _register_anndata(adata, data_registry_dict=data_registry)
    logger.info("Registered keys:{}".format(list(data_registry.keys())))
    _setup_summary_stats(adata, batch_key, labels_key, protein_expression_obsm_key)

    if copy:
        return adata


def assert_key_in_obs(adata, key):
    assert key in adata.obs.keys(), "{} is not a valid key for in adata.obs".format(key)


def _setup_labels(adata, labels_key):
    # checking labels
    if labels_key is None:
        logger.info("No label_key inputted, assuming all cells have same label")
        labels_key = "_scvi_labels"
        adata.obs[labels_key] = np.zeros(adata.shape[0], dtype=np.int64)
    else:
        assert_key_in_obs(adata, labels_key)
        logger.info('Using labels from adata.obs["{}"]'.format(labels_key))
        labels_key = _make_obs_column_categorical(
            adata, column_key=labels_key, alternate_column_key="_scvi_labels"
        )
    return labels_key


def _setup_batch(adata, batch_key):
    # checking batch
    # TODO Allow continuous batch information in the future
    if batch_key is None:
        logger.info("No batch_key inputted, assuming all cells are same batch")
        batch_key = "_scvi_batch"
        adata.obs[batch_key] = np.zeros(adata.shape[0], dtype=np.int64)
    else:
        assert_key_in_obs(adata, batch_key)
        logger.info('Using batches from adata.obs["{}"]'.format(batch_key))
        batch_key = _make_obs_column_categorical(
            adata, column_key=batch_key, alternate_column_key="_scvi_batch"
        )

    # make sure each batch contains enough cells
    unique, counts = np.unique(adata.obs[batch_key], return_counts=True)
    if np.min(counts) < 3:
        batch_val = unique[np.argmin(counts)]
        warnings.warn(
            "Batch {} has less than 3 cells. SCVI may not train properly.".format(
                batch_val
            )
        )
    # possible check for continuous?
    if len(unique) > (adata.shape[0] / 3):
        warnings.warn(
            "Is your batch continuous? SCVI doesn't support continuous batches yet."
        )
    return batch_key


def _make_obs_column_categorical(adata, column_key, alternate_column_key):
    """Makes the data in column_key in obs all categorical.
    if adata.obs[column_key] is not categorical, will categorize
    and save to .obs[alternate_column_key]
    """
    # check the datatype of data. if theyre not integers, make them ints
    user_data_dtype = adata.obs[column_key].dtype
    if user_data_dtype.name == "category":
        adata.obs[alternate_column_key] = (
            adata.obs[column_key].astype("category").cat.codes
        )
        column_key = alternate_column_key
    elif np.issubdtype(user_data_dtype, np.integer) is False:
        adata.obs[alternate_column_key] = (
            adata.obs[column_key].astype("category").cat.codes
        )
        column_key = alternate_column_key
    return column_key


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

    adata.uns["scvi_protein_names"] = protein_names

    # batch mask totalVI
    batch_mask = _get_batch_mask_protein_data(
        adata, protein_expression_obsm_key, batch_key
    )
    # check if it's actually needed
    if sum(sum([~b for b in batch_mask])) > 0:
        logger.info("Found batches with missing protein expression")
        adata.uns["totalvi_batch_mask"] = batch_mask
    return protein_expression_obsm_key


def _setup_X(adata, X_layers_key):
    # checking layers
    if X_layers_key is not None:
        assert (
            X_layers_key in adata.layers.keys()
        ), "{} is not a valid key in adata.layers".format(X_layers_key)
        logger.info('Using data from adata.layers["{}"]'.format(X_layers_key))
        X_loc = "layers"
        X_key = X_layers_key
        X = adata.layers[X_key]
    else:
        logger.info("Using data from adata.X")
        X_loc = "_X"
        X_key = None
        X = adata._X

    if _check_nonnegative_integers(X) is False:
        logger_data_loc = (
            "adata.X"
            if X_layers_key is None
            else "adata.layers[{}]".format(X_layers_key)
        )
        warnings.warn(
            "{} does not contain unnormalized count data. Are you sure this is what you want?".format(
                logger_data_loc
            )
        )
    if adata.shape[0] < adata.shape[1] < 1:
        warnings.warn("adata has more genes than cells. SCVI may not work properly.")

    return X_loc, X_key


def _setup_library_size(adata, batch_key, X_layers_key):
    # computes the library size per batch
    logger.info("Computing library size prior per batch")
    local_l_mean_key = "_scvi_local_l_mean"
    local_l_var_key = "_scvi_local_l_var"
    _compute_library_size_batch(
        adata,
        batch_key=batch_key,
        local_l_mean_key=local_l_mean_key,
        local_l_var_key=local_l_var_key,
        X_layers_key=X_layers_key,
    )
    return local_l_mean_key, local_l_var_key


def _setup_summary_stats(adata, batch_key, labels_key, protein_expression_obsm_key):
    n_batch = len(np.unique(adata.obs[batch_key]))
    n_cells = adata.shape[0]
    n_genes = adata.shape[1]
    n_labels = len(np.unique(adata.obs[labels_key]))

    summary_stats = {
        "n_batch": n_batch,
        "n_cells": n_cells,
        "n_genes": n_genes,
        "n_labels": n_labels,
    }
    if protein_expression_obsm_key is not None:
        summary_stats["n_proteins"] = adata.obsm[protein_expression_obsm_key].shape[1]
    adata.uns["scvi_summary_stats"] = summary_stats
    logger.info(
        "Successfully registered anndata object containing {} cells, {} genes, and {} batches.".format(
            n_cells, n_genes, n_batch
        )
    )
    return summary_stats
