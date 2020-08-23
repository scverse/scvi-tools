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
from pandas.api.types import CategoricalDtype

from scvi import _CONSTANTS

logger = logging.getLogger(__name__)


def get_from_registry(adata: anndata.AnnData, key: str) -> np.array:
    """Returns the object in AnnData associated with the key in ``adata.uns['_scvi']['data_registry']``

    Parameters
    ----------
    adata
        anndata object
    key
        key of object to get from ``adata.uns['_scvi]['data_registry']``

    Returns
    -------
    np.array containing the data requested

    Examples
    --------
    >>> import scvi
    >>> adata = scvi.dataset.cortex()
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
    if attr_key != "_use_raw":
        if isinstance(data, pd.DataFrame):
            data = data.loc[:, attr_key]
        else:
            data = data[attr_key]
    if isinstance(data, pd.Series):
        data = data.to_numpy().reshape(-1, 1)
    return data


def _transfer_categorical_mapping(adata1, adata2, scvi_key):
    """Makes the data in column_key in obs all categorical.
    if adata.obs[column_key] is not categorical, will categorize
    and save to .obs[alternate_column_key]

    Parameters:
    -----------
    adata
        anndata
    """
    categorical_mappings = adata1.uns["_scvi_categorical_mappings"]
    obs_key = categorical_mappings[scvi_key][0]
    obs_mapping = categorical_mappings[scvi_key][1]

    assert obs_key in adata2.obs.keys()

    # if theres more categories than previously registered
    assert len(obs_mapping) >= len(
        np.unique(adata2.obs[obs_key])
    ), "More batch categories than was registered"

    categorical_batch = adata2.obs[obs_key].astype(
        CategoricalDtype(categories=obs_mapping)
    )
    adata2.obs[scvi_key] = categorical_batch.cat.codes

    if np.min(adata2.obs[scvi_key]) == -1:
        logger.error(
            'Transfering category mapping for .obs["{}"] failed. It could be a datatype issue. Original adata category dtype: {}, New adata category dtype: {}'.format(
                obs_key, adata1.obs[obs_key][0].dtype, adata2.obs[obs_key][0].dtype
            )
        )
    if "_scvi_categorical_mappings" in adata2.uns.keys():
        adata2.uns["_scvi_categorical_mappings"][scvi_key] = [obs_key, obs_mapping]
    else:
        adata2.uns["_scvi_categorical_mappings"] = {
            scvi_key: [obs_key, obs_mapping],
        }


def transfer_anndata_setup(adata1, adata2, use_raw=False, X_layers_key=None):
    data_registry = adata1.uns["_scvi"]["data_registry"]

    # should we check that everything is the same size?
    if use_raw and X_layers_key:
        logging.warning(
            "use_raw and X_layers_key were both passed in. Defaulting to use_raw."
        )
    if use_raw:
        assert (
            get_from_registry(adata1, _CONSTANTS.X_KEY).shape[1]
            == adata2.raw.X.shape[1]
        )
    elif X_layers_key:
        assert (
            get_from_registry(adata1, _CONSTANTS.X_KEY).shape[1]
            == adata2.layers[X_layers_key].shape[1]
        )
    else:
        assert adata1.X.shape[1] == adata2.X.shape[1]

    has_protein = True if _CONSTANTS.PROTEIN_EXP_KEY in data_registry.keys() else False
    if has_protein is True:
        prev_protein_obsm_key = data_registry[_CONSTANTS.PROTEIN_EXP_KEY][1]
        if prev_protein_obsm_key not in adata2.obsm.keys():
            logger.warning(
                "Can't find {} in adata2.obsm for protein expressions. Ignoring protein expression.".format(
                    prev_protein_obsm_key
                )
            )
        else:
            assert (
                get_from_registry(adata1, _CONSTANTS.PROTEIN_EXP_KEY).shape[1]
                == adata2.obsm[prev_protein_obsm_key].shape[1]
            )
            protein_expression_obsm_key = prev_protein_obsm_key

    batch_key = "_scvi_batch"
    labels_key = "_scvi_labels"

    _transfer_categorical_mapping(adata1, adata2, batch_key)
    _transfer_categorical_mapping(adata1, adata2, labels_key)

    X_loc, X_key = _setup_X(adata2, X_layers_key, use_raw)
    local_l_mean_key, local_l_var_key = _setup_library_size(
        adata2, batch_key, X_layers_key, use_raw
    )
    data_registry = {
        _CONSTANTS.X_KEY: [X_loc, X_key],
        _CONSTANTS.BATCH_KEY: ["_obs", batch_key],
        _CONSTANTS.LOCAL_L_MEAN_KEY: ["_obs", local_l_mean_key],
        _CONSTANTS.LOCAL_L_VAR_KEY: ["_obs", local_l_var_key],
        _CONSTANTS.LABELS_KEY: ["_obs", labels_key],
    }

    if protein_expression_obsm_key is not None:
        protein_expression_obsm_key = _setup_protein_expression(
            adata2, protein_expression_obsm_key, None, batch_key
        )
        data_registry[_CONSTANTS.PROTEIN_EXP_KEY] = [
            "_obsm",
            protein_expression_obsm_key,
        ]
    # add the data_registry to anndata
    _register_anndata(adata2, data_registry_dict=data_registry)
    logger.info("Registered keys:{}".format(list(data_registry.keys())))
    _setup_summary_stats(adata2, batch_key, labels_key, protein_expression_obsm_key)


def setup_anndata(
    adata,
    batch_key: str = None,
    labels_key: str = None,
    use_raw: bool = False,
    X_layers_key: str = None,
    protein_expression_obsm_key: str = None,
    protein_names_uns_key: str = None,
    copy: bool = False,
) -> Optional[anndata.AnnData]:
    """Sets up AnnData object for scVI models.

    A mapping will be created between data fields used by scVI to their respective locations in adata
    This method will also compute the log mean and log variance per batch.

    None of the data in adata will be modified. Only adds fields to adata.

    Parameters
    ----------
    adata
        AnnData object containing raw counts
    batch_key
        key in adata.obs for batch information. If they're not integers, will automatically be converted into integer
        categories and saved to ``adata.obs['_scvi_batch']``. If ``None``, assigns the same batch to all the data.
    labels_key
        key in adata.obs for label information. If they're not integers, will automatically be converted into integer
        categories and saved to ``adata.obs['_scvi_labels']``. If ``None``, assigns the same label to all the data.
    X_layers_key
        if not None, uses this as the key in ``adata.layers`` for raw count data.
    protein_expression_obsm_key
        key in ``adata.obsm`` for protein expression data, Required for TotalVI.
    protein_names_uns_key
        key in ``adata.uns`` for protein names. If None, will use the column names of ``adata.obsm[protein_expression_obsm_key]``
        if it is a pandas dataframe, else will assign sequential names to proteins. Only relavent but not required for TotalVI.
    copy
        if True, a copy of adata is returned.

    Returns
    -------
    if ``copy``,  will return ``anndata.AnnData`` else adds the following fields to adata:

    .uns['scvi_data_registy']
        dictionary mapping data fields used by scVI to their respective locations in adata
    .uns['scvi_summary_stats']
        dictionary of summary statistics for adata
    .obs['_local_l_mean']
        per batch library size mean
    .obs['_local_l_var']
        per batch library size variance

    if no ``batch_key`` or ``labels_key`` was provided, or if they were not encoded as integers, will also add the following:

    adata.obs['_scvi_labels']
        labels encoded as integers
    adata.obs['_scvi_batch']
        batch encoded as integers

    Examples
    --------

    Example setting up a scanpy dataset with random gene data and no batch nor label information

    >>> import scanpy
    >>> import scvi
    >>> import numpy
    >>> adata = scanpy.datasets.blobs()
    >>> # filter cells and run preprocessing before setup_anndata
    >>> scanpy.pp.filter_cells(adata, min_counts = 0)
    >>> # since no batch_key nor labels_key was passed, setup_anndata() will assume all cells have the same batch and label
    >>> setup_anndata(adata)
    [2020-07-31 12:54:15,293] INFO - scvi.dataset._anndata | No batch_key inputted, assuming all cells are same batch
    [2020-07-31 12:54:15,296] INFO - scvi.dataset._anndata | No label_key inputted, assuming all cells have same label
    [2020-07-31 12:54:15,298] INFO - scvi.dataset._anndata | Using data from adata.X
    [2020-07-31 12:54:15,300] INFO - scvi.dataset._anndata | Computing library size prior per batch
    [2020-07-31 12:54:15,309] INFO - scvi.dataset._anndata | Registered keys:['X', 'batch_indices', 'local_l_mean', 'local_l_var', 'labels']
    [2020-07-31 12:54:15,312] INFO - scvi.dataset._anndata | Successfully registered anndata object containing 421 cells, 11 genes, and 1 batches.
    >>> # see registered scVI fields and their respective locations in adata
    >>> adata.uns['_scvi']['scvi_data_registry']
    {'X': ['_X', None],
    'batch_indices': ['obs', '_scvi_batch'],
    'local_l_mean': ['obs', '_scvi_local_l_mean'],
    'local_l_var': ['obs', '_scvi_local_l_var'],
    'labels': ['obs', '_scvi_labels']}
    >>> # summary statistics can be used to spot check setup
    >>> adata.uns['scvi_summary_stats']
    {'n_batch': 1, 'n_cells': 421, 'n_genes': 11, 'n_labels': 1}

    Example setting up scanpy dataset with random gene data, batch, and protein expression

    >>> import scanpy
    >>> import scvi
    >>> import numpy
    >>> adata = scanpy.datasets.blobs()
    >>> # filter cells and run preprocessing before setup_anndata
    >>> scanpy.pp.filter_cells(adata, min_counts = 0)
    >>> adata
    AnnData object with n_obs × n_vars = 421 × 11
        obs: 'blobs'
    >>> # generate random protein expression data
    >>> adata.obsm['rand_protein_exp'] = numpy.random.randint(0, 100, (adata.n_obs, 20))
    >>> # setup adata with batch info from adata.obs['blobs'] and protein expressions from adata.obsm['rand_protein_exp']
    >>> scvi.dataset.setup_anndata(adata, batch_key='blobs', protein_expression_obsm_key='rand_protein_exp' )
    [2020-07-31 12:44:30,421] INFO - scvi.dataset._anndata | Using batches from adata.obs["blobs"]
    [2020-07-31 12:44:30,425] INFO - scvi.dataset._anndata | No label_key inputted, assuming all cells have same label
    [2020-07-31 12:44:30,427] INFO - scvi.dataset._anndata | Using data from adata.X
    [2020-07-31 12:44:30,428] INFO - scvi.dataset._anndata | Computing library size prior per batch
    [2020-07-31 12:44:30,442] INFO - scvi.dataset._anndata | Using protein expression from adata.obsm['rand_protein_exp']
    [2020-07-31 12:44:30,443] INFO - scvi.dataset._anndata | Generating sequential protein names
    [2020-07-31 12:44:30,445] INFO - scvi.dataset._anndata | Registered keys:['X', 'batch_indices', 'local_l_mean', 'local_l_var', 'labels', 'protein_expression']
    [2020-07-31 12:44:30,446] INFO - scvi.dataset._anndata | Successfully registered anndata object containing 421 cells, 11 genes, and 4 batches.
    >>> # see registered scVI fields and their respective locations in adata
    >>> adata.uns['_scvi]['data_registry']
    {'X': ['_X', None],
    'batch_indices': ['obs', '_scvi_batch'],
    'local_l_mean': ['obs', '_scvi_local_l_mean'],
    'local_l_var': ['obs', '_scvi_local_l_var'],
    'labels': ['obs', '_scvi_labels'],
    'protein_expression': ['obsm', 'rand_protein_exp']}
    >>> # summary statistics can be used to spot check setup
    >>> adata.uns['scvi_summary_stats']
    {'n_batch': 4, 'n_cells': 421, 'n_genes': 11, 'n_labels': 1, 'n_proteins': 20}

    Example with counts in a layer in AnnData and protein expression data as a pandas DataFrame

    >>> import scvi
    >>> # load a dataset
    >>> adata = scvi.dataset.pbmcs_10x_cite_seq(run_setup_anndata = False)
    >>> adata.layers['raw_counts'] = adata.X.copy()
    >>> # keys in adata.obs
    >>> adata.obs.keys()
    Index(['n_genes', 'percent_mito', 'n_counts', 'batch', 'labels'], dtype='object')
    >>> # adata.obsm['protein_expression'] contains expression data as well as protein names
    >>> adata.obsm['protein_expression'].head()
                            CD3_TotalSeqB  CD4_TotalSeqB  ...  TIGIT_TotalSeqB  CD127_TotalSeqB
    index                                               ...
    AAACCCAAGATTGTGA-1-0             18            138  ...                4                7
    AAACCCACATCGGTTA-1-0             30            119  ...                9                8
    AAACCCAGTACCGCGT-1-0             18            207  ...               11               12
    AAACCCAGTATCGAAA-1-0             18             11  ...               59               16
    AAACCCAGTCGTCATA-1-0              5             14  ...               76               17
    [5 rows x 14 columns]
    >>> # setup anndata with batch and labels information from adata.obs and protein expressions from adata.obsm
    >>> setup_anndata(adata, batch_key='batch', labels_key='labels', protein_expression_obsm_key='protein_expression', X_layers_key='raw_counts')
    [2020-07-31 13:04:28,129] INFO - scvi.dataset._anndata | Using batches from adata.obs["batch"]
    [2020-07-31 13:04:28,135] INFO - scvi.dataset._anndata | Using labels from adata.obs["labels"]
    [2020-07-31 13:04:28,136] INFO - scvi.dataset._anndata | Using data from adata.layers["raw_counts"]
    [2020-07-31 13:04:30,018] INFO - scvi.dataset._anndata | Computing library size prior per batch
    [2020-07-31 13:04:30,664] INFO - scvi.dataset._anndata | Using protein expression from adata.obsm['protein_expression']
    [2020-07-31 13:04:30,669] INFO - scvi.dataset._anndata | Using protein names from columns of adata.obsm['protein_expression']
    [2020-07-31 13:04:30,672] INFO - scvi.dataset._anndata | Registered keys:['X', 'batch_indices', 'local_l_mean', 'local_l_var', 'labels', 'protein_expression']
    [2020-07-31 13:04:30,673] INFO - scvi.dataset._anndata | Successfully registered anndata object containing 10849 cells, 15792 genes, and 2 batches.
    >>> # see registered scVI fields and their respective loaations in adata
    >>> adata.uns['scvi']['data_registry']
    {'X': ['layers', 'raw_counts'],
    'batch_indices': ['obs', 'batch'],
    'local_l_mean': ['obs', '_scvi_local_l_mean'],
    'local_l_var': ['obs', '_scvi_local_l_var'],
    'labels': ['obs', 'labels'],
    'protein_expression': ['obsm', 'protein_expression']}
    >>> # summary statistics can be used to spot check setup
    >>> adata.uns['scvi_summary_stats']
    {'n_batch': 2,
    'n_cells': 10849,
    'n_genes': 15792,
    'n_labels': 1,
    'n_proteins': 14}
    >>> # notice that this is the same as the columns of adata.obsm['protein_expression']
    >>> adata.uns['scvi_protein_names']
    ['CD3_TotalSeqB',
    'CD4_TotalSeqB',
    'CD8a_TotalSeqB',
    'CD14_TotalSeqB',
    'CD15_TotalSeqB',
    'CD16_TotalSeqB',
    'CD56_TotalSeqB',
    'CD19_TotalSeqB',
    'CD25_TotalSeqB',
    'CD45RA_TotalSeqB',
    'CD45RO_TotalSeqB',
    'PD-1_TotalSeqB',
    'TIGIT_TotalSeqB',
    'CD127_TotalSeqB']

    """
    if adata.is_view:
        raise ValueError("adata cannot be a view of an AnnData object.")

    if copy:
        adata = adata.copy()
    adata.uns["_scvi"] = {}

    batch_key = _setup_batch(adata, batch_key)
    labels_key = _setup_labels(adata, labels_key)
    X_loc, X_key = _setup_X(adata, X_layers_key, use_raw)
    local_l_mean_key, local_l_var_key = _setup_library_size(
        adata, batch_key, X_layers_key, use_raw
    )

    data_registry = {
        _CONSTANTS.X_KEY: {"attr_name": X_loc, "attr_key": X_key},
        _CONSTANTS.BATCH_KEY: {"attr_name": "_obs", "attr_key": batch_key},
        _CONSTANTS.LOCAL_L_MEAN_KEY: {
            "attr_name": "_obs",
            "attr_key": local_l_mean_key,
        },
        _CONSTANTS.LOCAL_L_VAR_KEY: {"attr_name": "_obs", "attr_key": local_l_var_key},
        _CONSTANTS.LABELS_KEY: {"attr_name": "_obs", "attr_key": labels_key},
    }

    if protein_expression_obsm_key is not None:
        protein_expression_obsm_key = _setup_protein_expression(
            adata, protein_expression_obsm_key, protein_names_uns_key, batch_key
        )
        data_registry[_CONSTANTS.PROTEIN_EXP_KEY] = {
            "attr_name": "_obsm",
            "attr_key": protein_expression_obsm_key,
        }
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
        alt_key = "_scvi_{}".format(labels_key)
        labels_key = _make_obs_column_categorical(
            adata, column_key=labels_key, alternate_column_key=alt_key
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
        alt_key = "_scvi_{}".format(batch_key)
        batch_key = _make_obs_column_categorical(
            adata, column_key=batch_key, alternate_column_key=alt_key
        )
    return batch_key


def _make_obs_column_categorical(adata, column_key, alternate_column_key):
    """Makes the data in column_key in obs all categorical.
    if adata.obs[column_key] is not categorical, will categorize
    and save to .obs[alternate_column_key]
    """
    categorical_obs = adata.obs[column_key].astype("category")
    adata.obs[alternate_column_key] = categorical_obs.cat.codes
    mapping = categorical_obs.cat.categories
    if "categorical_mappings" in adata.uns["_scvi"].keys():
        adata.uns["_scvi"]["categorical_mappings"][alternate_column_key] = [
            column_key,
            mapping,
        ]
    else:
        adata.uns["_scvi"]["categorical_mappings"] = {
            alternate_column_key: [column_key, mapping]
        }
    column_key = alternate_column_key

    # make sure each category contains enough cells
    unique, counts = np.unique(adata.obs[column_key], return_counts=True)
    if np.min(counts) < 3:
        category = unique[np.argmin(counts)]
        warnings.warn(
            "Category {} in adata.obs['{}'] has fewer than 3 cells. SCVI may not train properly.".format(
                category, column_key
            )
        )
    # possible check for continuous?
    if len(unique) > (adata.shape[0] / 3):
        warnings.warn(
            "Is adata.obs['{}'] continuous? SCVI doesn't support continuous obs yet."
        )
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
    if np.sum([~b for b in batch_mask]) > 0:
        logger.info("Found batches with missing protein expression")
        adata.uns["totalvi_batch_mask"] = batch_mask
    return protein_expression_obsm_key


def _setup_X(adata, X_layers_key, use_raw):
    if use_raw and X_layers_key:
        logging.warning(
            "use_raw and X_layers_key were both passed in. Defaulting to use_raw."
        )

    # checking layers
    if use_raw:
        assert adata.raw is not None, "use_raw is True but adata.raw is None"
        logger.info("Using data from adata.raw.X")
        X_loc = "raw"
        X_key = "X"
        X = adata.raw._X
    elif X_layers_key is not None:
        assert (
            X_layers_key in adata.layers.keys()
        ), "{} is not a valid key in adata.layers".format(X_layers_key)
        logger.info('Using data from adata.layers["{}"]'.format(X_layers_key))
        X_loc = "_layers"
        X_key = X_layers_key
        X = adata.layers[X_key]
    else:
        logger.info("Using data from adata.X")
        X_loc = "_X"
        X_key = "_use_raw"
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
    if adata.shape[0] / adata.shape[1] < 1:
        warnings.warn("adata has more genes than cells. SCVI may not work properly.")

    return X_loc, X_key


def _setup_library_size(adata, batch_key, X_layers_key, use_raw):
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
        use_raw=use_raw,
    )
    return local_l_mean_key, local_l_var_key


def _setup_summary_stats(adata, batch_key, labels_key, protein_expression_obsm_key):
    n_batch = len(np.unique(adata.obs[batch_key]))
    n_cells = adata.shape[0]
    n_genes = adata.shape[1]
    n_labels = len(np.unique(adata.obs[labels_key]))

    if protein_expression_obsm_key is not None:
        n_proteins = adata.obsm[protein_expression_obsm_key].shape[1]
    else:
        n_proteins = 0

    summary_stats = {
        "n_batch": n_batch,
        "n_cells": n_cells,
        "n_genes": n_genes,
        "n_labels": n_labels,
        "n_proteins": n_proteins,
    }
    adata.uns["_scvi"]["summary_stats"] = summary_stats
    logger.info(
        "Successfully registered anndata object containing {} cells, {} genes, {} batches, and {} proteins.".format(
            n_cells, n_genes, n_batch, n_proteins
        )
    )
    return summary_stats


def _register_anndata(adata, data_registry_dict: Dict[str, Tuple[str, str]]):
    """Registers the AnnData object by adding data_registry_dict to adata.uns['_scvi']['data_registry']

    Format of data_registry_dict is: {<scvi_key>: (<anndata dataframe>, <dataframe key> )}

    Parameters
    ----------
    adata
        anndata object
    data_registry_dict
        dictionary mapping keys used by scvi models to their respective location in adata.

    Examples
    --------
    >>> data_dict = {"batch" :("obs", "batch_idx"), "X": ("_X", None)}
    >>> _register_anndata(adata, data_dict)

    """
    # check doesnt work for adata.raw since adata.raw is not a dict

    # for df, df_key in data_registry_dict.values():
    #     if df_key is not None:
    #         assert df_key in getattr(
    #             adata, df
    #         ), "anndata.{} has no attribute '{}'".format(df, df_key)
    #     else:
    #         assert hasattr(adata, df) is True, "anndata has no attribute '{}'".format(
    #             df_key
    #         )
    adata.uns["_scvi"]["data_registry"] = copy.copy(data_registry_dict)
