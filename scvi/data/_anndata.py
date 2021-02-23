import logging
import os
import pickle
import sys
import warnings
from typing import List, Optional, Union

import anndata
import numpy as np
import pandas as pd
import rich
from anndata._core.anndata import AnnData
from pandas.api.types import CategoricalDtype
from rich.console import Console

import scvi
from scvi import _CONSTANTS
from scvi._compat import Literal

from ._anndatarecorder import AnnDataRecorder
from ._utils import (
    _check_nonnegative_integers,
    _compute_library_size_batch,
    _get_batch_mask_protein_data,
)

logger = logging.getLogger(__name__)


def setup_anndata(
    adata: anndata.AnnData,
    batch_key: Optional[str] = None,
    labels_key: Optional[str] = None,
    layer: Optional[str] = None,
    protein_expression_obsm_key: Optional[str] = None,
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

    recorder = AnnDataRecorder(adata)
    recorder.setup_batch(batch_key)
    recorder.setup_labels(labels_key)
    recorder.setup_x(layer)

    local_l_mean_key, local_l_var_key = _setup_library_size(adata, batch_key, layer)

    # TODO: move this code to be model specific
    data_registry = adata.uns["_scvi"]["data_registry"]
    data_registry.update(
        {
            _CONSTANTS.LOCAL_L_MEAN_KEY: {
                "attr_name": "obs",
                "attr_key": local_l_mean_key,
            },
            _CONSTANTS.LOCAL_L_VAR_KEY: {
                "attr_name": "obs",
                "attr_key": local_l_var_key,
            },
        }
    )

    if protein_expression_obsm_key is not None:
        protein_expression_obsm_key = _setup_protein_expression(
            recorder, protein_expression_obsm_key, batch_key
        )
    else:
        recorder.add_to_summary_stats("n_proteins", 0)

    if categorical_covariate_keys is not None:
        recorder.setup_categorical_obs_key_iterable(
            categorical_covariate_keys=categorical_covariate_keys,
            recorded_key="_scvi_extra_categoricals",
            registry_key=_CONSTANTS.CAT_COVS_KEY,
        )
        n_cat_keys = len(categorical_covariate_keys)
    else:
        n_cat_keys = 0
    recorder.add_to_summary_stats("n_categorical_covs", n_cat_keys)

    if continuous_covariate_keys is not None:
        recorder.setup_continuous_obs_key_iterable(
            continuous_covariate_keys=continuous_covariate_keys,
            recorded_key="_scvi_extra_continuous",
            registry_key=_CONSTANTS.CONT_COVS_KEY,
        )
        n_cont_keys = len(continuous_covariate_keys)
    else:
        n_cont_keys = 0
    recorder.add_to_summary_stats("n_continuous_covs", n_cont_keys)

    recorder.set_setup_dict()

    ss = recorder.setup_dict["summary_stats"]
    logger.info(
        "Successfully registered anndata object containing {} cells, {} vars, "
        "{} batches, {} labels, and {} proteins. Also registered {} extra categorical "
        "covariates and {} extra continuous covariates.".format(
            ss["n_cells"],
            ss["n_vars"],
            ss["n_batch"],
            ss["n_labels"],
            ss["n_proteins"],
            ss["n_categorical_covs"],
            ss["n_continuous_covs"],
        )
    )
    logger.info("Please do not further modify adata until model is trained.")

    if copy:
        return adata


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
    target_recorder = AnnDataRecorder(adata_target)
    target_setup_dict = target_recorder.setup_dict

    if isinstance(adata_source, anndata.AnnData):
        _scvi_dict = adata_source.uns["_scvi"]
    else:
        _scvi_dict = adata_source
    data_registry = _scvi_dict["data_registry"]
    summary_stats = _scvi_dict["summary_stats"]

    target_recorder.setup_dict["data_registry"] = _scvi_dict["data_registry"].copy()

    # transfer version
    target_setup_dict["scvi_version"] = _scvi_dict["scvi_version"]

    # transfer X
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
    target_recorder.setup_x(layer)

    # transfer categorical obs keys
    categorical_mappings = _scvi_dict["categorical_mappings"]
    _transfer_categorical_obs(target_recorder, categorical_mappings, extend_categories)

    batch_key = "_scvi_batch"
    labels_key = "_scvi_labels"

    # library size
    local_l_mean_key, local_l_var_key = _setup_library_size(
        adata_target, batch_key, layer
    )

    # transfer protein_expression
    protein_expression_obsm_key = _transfer_protein_expression(
        _scvi_dict, target_recorder, batch_key
    )

    # transfer extra categorical covs
    # TODO: loop over every key in _scvi_dict["categorical_obsm_keys"]
    has_cat_cov = True if _CONSTANTS.CAT_COVS_KEY in data_registry.keys() else False
    if has_cat_cov:
        source_cat_dict = _scvi_dict["categorical_obsm_keys"][
            "_scvi_extra_categoricals"
        ]["mappings"].copy()
        # extend categories
        if extend_categories:
            for key, mapping in source_cat_dict:
                for c in np.unique(adata_target.obs[key]):
                    if c not in mapping:
                        mapping = np.concatenate([mapping, [c]])
                source_cat_dict[key] = mapping
        # use ["keys"] to maintain correct order
        target_recorder.setup_categorical_obs_key_iterable(
            categorical_covariate_keys=_scvi_dict["categorical_obsm_keys"][
                "_scvi_extra_categoricals"
            ]["keys"],
            recorded_key="_scvi_extra_categoricals",
            registry_key=_CONSTANTS.CAT_COVS_KEY,
            category_dict=source_cat_dict,
        )
    else:
        source_cat_dict = None

    # transfer extra continuous covs
    # TODO: loop over every key in _scvi_dict["categorical_obsm_keys"]
    has_cont_cov = True if _CONSTANTS.CONT_COVS_KEY in data_registry.keys() else False
    if has_cont_cov:
        obs_keys_names = _scvi_dict["continuous_obsm_keys"]["_scvi_extra_continuous"][
            "keys"
        ]
        target_recorder.setup_continuous_obs_key_iterable(
            continuous_covariate_keys=list(obs_keys_names),
            recorded_key="_scvi_extra_continuous",
            registry_key=_CONSTANTS.CONT_COVS_KEY,
        )
    else:
        obs_keys_names = None

    target_recorder.set_setup_dict()

    _setup_summary_stats(
        adata_target,
        batch_key,
        labels_key,
        protein_expression_obsm_key,
        source_cat_dict,
        obs_keys_names,
    )


def _transfer_categorical_obs(target_recorder, categorical_mappings, extend_categories):
    adata_target = target_recorder.adata
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
        target_recorder._make_obs_column_categorical(original_key, key, cat_dtype)


def _transfer_protein_expression(_scvi_dict, target_recorder, batch_key):
    adata_target = target_recorder.adata
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

            _setup_protein_expression(
                target_recorder, protein_expression_obsm_key, batch_key
            )
    else:
        protein_expression_obsm_key = None


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


def _setup_protein_expression(recorder, protein_expression_obsm_key, batch_key):
    adata = recorder.adata
    pro_exp = adata.obsm[protein_expression_obsm_key]
    if _check_nonnegative_integers(pro_exp) is False:
        warnings.warn(
            "adata.obsm[{}] does not contain unnormalized count data. Are you sure this is what you want?".format(
                protein_expression_obsm_key
            )
        )
    recorder.setup_obsm_key(protein_expression_obsm_key, _CONSTANTS.PROTEIN_EXP_KEY)
    logger.info(
        "Using protein expression from adata.obsm['{}']".format(
            protein_expression_obsm_key
        )
    )
    # setup protein names
    if isinstance(adata.obsm[protein_expression_obsm_key], pd.DataFrame):
        logger.info(
            "Using protein names from columns of adata.obsm['{}']".format(
                protein_expression_obsm_key
            )
        )
        protein_names = list(adata.obsm[protein_expression_obsm_key].columns)
    else:
        logger.info("Generating sequential protein names")
        protein_names = np.arange(adata.obsm[protein_expression_obsm_key].shape[1])

    recorder.setup_dict["protein_names"] = protein_names
    recorder.add_to_summary_stats("n_proteins", len(protein_names))

    # batch mask totalVI
    batch_mask = _get_batch_mask_protein_data(
        adata, protein_expression_obsm_key, batch_key
    )

    # check if it's actually needed
    if np.sum([~b[1] for b in batch_mask.items()]) > 0:
        logger.info("Found batches with missing protein expression")
        recorder.setup_dict["protein_names"]["totalvi_batch_mask"] = batch_mask
    return protein_expression_obsm_key


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


def _check_anndata_setup_equivalence(
    adata_source: Union[AnnData, dict], adata_target: AnnData
) -> bool:
    """
    Checks if target setup is equivalent to source.

    Parameters
    ----------
    adata_source
        Either AnnData already setup or scvi_setup_dict as the source
    adata_target
        Target AnnData to check setup equivalence

    Returns
    -------
    Whether the adata_target should be run through `transfer_anndata_setup`
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

    return transfer_setup


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
