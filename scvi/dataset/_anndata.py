import copy
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
from scvi.dataset._constants import (
    _X_KEY,
    _BATCH_KEY,
    _LOCAL_L_MEAN_KEY,
    _LOCAL_L_VAR_KEY,
    _LABELS_KEY,
    _PROTEIN_EXP_KEY,
)

logger = logging.getLogger(__name__)


def _register_anndata(adata, data_registry_dict: Dict[str, Tuple[str, str]]):
    """Registers the AnnData object by adding data_registry_dict to adata.uns

    Format is: {<scvi_key>: (<anndata dataframe>, <dataframe key> )}
    Example:
    {"batch" :("obs", "batch_idx")}
    {"X": (None, "X")}

    Parameters
    ----------
    adata
        anndata object
    data_registry_dict
        dictionary mapping keys used by scvi models to their respective location in adata.

    """
    for df, df_key in data_registry_dict.values():
        if df is not None:
            assert df_key in getattr(
                adata, df
            ), "anndata.{} has no attribute '{}'".format(df, df_key)
        else:
            assert (
                hasattr(adata, df_key) is True
            ), "anndata has no attribute '{}'".format(df_key)
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
    # assert "scvi_data_registry" in adata.uns.keys(), "AnnData was never registered"
    data_loc = adata.uns["scvi_data_registry"][key]
    df, df_key = data_loc[0], data_loc[1]
    if df == "":
        df = None
    data = getattr(adata, df)[df_key] if df is not None else getattr(adata, df_key)
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

    # if dataset takes less than 1gb data, can make it dense
    if copy:
        adata = adata.copy()

    # checking layers
    if X_layers_key is not None and X_layers_key not in adata.layers.keys():
        raise ValueError("{} is not a valid key in adata.layers".format(X_layers_key))

    X = adata.layers[X_layers_key] if X_layers_key is not None else adata.X
    if _check_nonnegative_integers(X) is False:
        logger_data_loc = (
            "adata.X"
            if X_layers_key is None
            else "adata.layers[{}]".format(X_layers_key)
        )
        logger.warning(
            "{} does not contain unnormalized count data. Are you sure this is what you want?".format(
                logger_data_loc
            )
        )
    if X_layers_key is not None:
        logger.info('Using data from adata.layers["{}"]'.format(X_layers_key))
    else:
        logger.info("Using data from adata.X")

    # checking batch
    if batch_key is None:
        logger.info("No batch_key inputted, assuming all cells are same batch")
        batch_key = "_scvi_batch"
        adata.obs[batch_key] = np.zeros(adata.shape[0])
    else:
        assert (
            batch_key in adata.obs.keys()
        ), "{} is not a valid key in adata.obs".format(batch_key)
        logger.info('Using batches from adata.obs["{}"]'.format(batch_key))

    # check the datatype of batches. if theyre not integers, make them ints
    user_batch_dtype = adata.obs[batch_key].dtype
    if user_batch_dtype.name == "category":
        adata.obs["_scvi_batch"] = adata.obs[batch_key].astype("category").cat.codes
        batch_key = "_scvi_batch"
    elif np.issubdtype(user_batch_dtype, np.integer) is False:
        adata.obs["_scvi_batch"] = adata.obs[batch_key].astype("category").cat.codes
        batch_key = "_scvi_batch"

    if labels_key is None:
        logger.info("No label_key inputted, assuming all cells have same label")
        labels_key = "_scvi_labels"
        adata.obs[labels_key] = np.zeros(adata.shape[0])
    else:
        assert (
            labels_key in adata.obs.keys()
        ), "{} is not a valid key for in adata.obs".format(labels_key)
        logger.info('Using labels from adata.obs["{}"]'.format(labels_key))

    # check the datatype of labels. if theyre not integers, make them ints
    user_labels_dtype = adata.obs[labels_key].dtype
    if user_labels_dtype.name == "category":
        adata.obs["_scvi_labels"] = adata.obs[labels_key].astype("category").cat.codes
        labels_key = "_scvi_labels"
    elif np.issubdtype(user_labels_dtype, np.integer) is False:
        adata.obs["_scvi_labels"] = adata.obs[labels_key].astype("category").cat.codes
        labels_key = "_scvi_labels"

    # computes the library size per batch
    local_l_mean_key = "_scvi_local_l_mean"
    local_l_var_key = "_scvi_local_l_var"

    logger.info("Computing library size prior per batch")

    _compute_library_size_batch(
        adata,
        batch_key=batch_key,
        local_l_mean_key=local_l_mean_key,
        local_l_var_key=local_l_var_key,
        X_layers_key=X_layers_key,
    )

    if X_layers_key is None:
        X_loc = None
        X_key = "_X"
    else:
        X_loc = "layers"
        X_key = X_layers_key

    data_registry = {
        _X_KEY: [X_loc, X_key],
        _BATCH_KEY: ["obs", batch_key],
        _LOCAL_L_MEAN_KEY: ["obs", local_l_mean_key],
        _LOCAL_L_VAR_KEY: ["obs", local_l_var_key],
        _LABELS_KEY: ["obs", labels_key],
    }

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
        assert (
            protein_expression_obsm_key in adata.obsm.keys()
        ), "{} is not a valid key in adata.obsm".format(protein_expression_obsm_key)
        pro_exp = adata.obsm[protein_expression_obsm_key]
        if _check_nonnegative_integers(pro_exp) is False:
            logger.warning(
                "adata.obsm[{}] does not contain unnormalized count data. Are you sure this is what you want?".format(
                    protein_expression_obsm_key
                )
            )
        logger.info(
            "Using protein expression from adata.obsm['{}']".format(
                protein_expression_obsm_key
            )
        )

        data_registry[_PROTEIN_EXP_KEY] = ["obsm", protein_expression_obsm_key]
        summary_stats["n_proteins"] = adata.obsm[protein_expression_obsm_key].shape[1]
        if protein_names_uns_key is None and isinstance(
            adata.obsm[protein_expression_obsm_key], pd.DataFrame
        ):
            summary_stats["protein_names"] = list(
                adata.obsm[protein_expression_obsm_key].columns
            )
            logger.info(
                "Using protein names from columns of adata.obsm['{}']".format(
                    protein_expression_obsm_key
                )
            )
        elif protein_names_uns_key is not None:
            summary_stats["protein_names"] = adata.uns[protein_names_uns_key]
            logger.info(
                "Using protein names from adata.uns['{}']".format(protein_names_uns_key)
            )
        else:
            summary_stats["protein_names"] = np.arange(
                adata.obsm[protein_expression_obsm_key].shape[1]
            )
            logger.info("Generating sequential protein names")
        # batch mask totalVI
        batch_mask = _get_batch_mask_protein_data(
            adata, protein_expression_obsm_key, batch_key
        )
        # check if it's actually needed
        if sum(sum([~b for b in batch_mask])) > 0:
            logger.info("Found batches with missing protein expression")
            adata.uns["totalvi_batch_mask"] = batch_mask

    _register_anndata(adata, data_registry_dict=data_registry)

    logger.info(
        "Successfully registered anndata object containing {} cells, {} genes, and {} batches \nRegistered keys:{}".format(
            n_cells, n_genes, n_batch, list(data_registry.keys())
        )
    )
    adata.uns["scvi_summary_stats"] = summary_stats
    if copy:
        return adata
