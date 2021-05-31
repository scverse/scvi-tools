import logging
import os
import pickle
import warnings
from collections.abc import Iterable as IterableClass
from typing import List, Optional, Union

import numpy as np
import pandas as pd
import torch
from anndata import AnnData, read

from scvi._compat import Literal
from scvi.utils import DifferentialComputation, track

logger = logging.getLogger(__name__)


def _load_saved_files(
    dir_path: str,
    load_adata: bool,
    map_location: Optional[Literal["cpu", "cuda"]] = None,
):
    """Helper to load saved files."""
    setup_dict_path = os.path.join(dir_path, "attr.pkl")
    adata_path = os.path.join(dir_path, "adata.h5ad")
    varnames_path = os.path.join(dir_path, "var_names.csv")
    model_path = os.path.join(dir_path, "model_params.pt")

    if os.path.exists(adata_path) and load_adata:
        adata = read(adata_path)
    elif not os.path.exists(adata_path) and load_adata:
        raise ValueError("Save path contains no saved anndata and no adata was passed.")
    else:
        adata = None

    var_names = np.genfromtxt(varnames_path, delimiter=",", dtype=str)

    with open(setup_dict_path, "rb") as handle:
        attr_dict = pickle.load(handle)

    scvi_setup_dict = attr_dict.pop("scvi_setup_dict_")

    model_state_dict = torch.load(model_path, map_location=map_location)

    return scvi_setup_dict, attr_dict, var_names, model_state_dict, adata


def _initialize_model(cls, adata, attr_dict):
    """Helper to initialize a model."""
    if "init_params_" not in attr_dict.keys():
        raise ValueError(
            "No init_params_ were saved by the model. Check out the "
            "developers guide if creating custom models."
        )
    # get the parameters for the class init signiture
    init_params = attr_dict.pop("init_params_")

    # new saving and loading, enable backwards compatibility
    if "non_kwargs" in init_params.keys():
        # grab all the parameters execept for kwargs (is a dict)
        non_kwargs = init_params["non_kwargs"]
        kwargs = init_params["kwargs"]

        # expand out kwargs
        kwargs = {k: v for (i, j) in kwargs.items() for (k, v) in j.items()}
    else:
        # grab all the parameters execept for kwargs (is a dict)
        non_kwargs = {k: v for k, v in init_params.items() if not isinstance(v, dict)}
        kwargs = {k: v for k, v in init_params.items() if isinstance(v, dict)}
        kwargs = {k: v for (i, j) in kwargs.items() for (k, v) in j.items()}
        non_kwargs.pop("use_cuda")

    model = cls(adata, **non_kwargs, **kwargs)
    return model


def _validate_var_names(adata, source_var_names):

    user_var_names = adata.var_names.astype(str)
    if not np.array_equal(source_var_names, user_var_names):
        warnings.warn(
            "var_names for adata passed in does not match var_names of "
            "adata used to train the model. For valid results, the vars "
            "need to be the same and in the same order as the adata used to train the model."
        )


def _prepare_obs(
    idx1: Union[List[bool], np.ndarray, str],
    idx2: Union[List[bool], np.ndarray, str],
    adata: AnnData,
):
    """
    Construct an array used for masking.

    Given population identifiers `idx1` and potentially `idx2`,
    this function creates an array `obs_col` that identifies both populations
    for observations contained in `adata`.
    In particular, `obs_col` will take values `group1` (resp. `group2`)
    for `idx1` (resp `idx2`).

    Parameters
    ----------
    idx1
        Can be of three types. First, it can corresponds to a boolean mask that
        has the same shape as adata. It can also corresponds to a list of indices.
        Last, it can correspond to string query of adata.obs columns.
    idx2
        Same as above
    adata
        Anndata
    """

    def ravel_idx(my_idx, obs_df):
        return (
            obs_df.index.isin(obs_df.query(my_idx).index)
            if isinstance(my_idx, str)
            else np.asarray(my_idx).ravel()
        )

    obs_df = adata.obs
    idx1 = ravel_idx(idx1, obs_df)
    g1_key = "one"
    obs_col = np.array(["None"] * adata.shape[0], dtype=str)
    obs_col[idx1] = g1_key
    group1 = [g1_key]
    group2 = None if idx2 is None else "two"
    if idx2 is not None:
        idx2 = ravel_idx(idx2, obs_df)
        obs_col[idx2] = group2
    if (obs_col[idx1].shape[0] == 0) or (obs_col[idx2].shape[0] == 0):
        raise ValueError("One of idx1 or idx2 has size zero.")
    return obs_col, group1, group2


def _de_core(
    adata,
    model_fn,
    groupby,
    group1,
    group2,
    idx1,
    idx2,
    all_stats,
    all_stats_fn,
    col_names,
    mode,
    batchid1,
    batchid2,
    delta,
    batch_correction,
    fdr,
    silent,
    **kwargs
):
    """Internal function for DE interface."""
    if group1 is None and idx1 is None:
        group1 = adata.obs[groupby].astype("category").cat.categories.tolist()
        if len(group1) == 1:
            raise ValueError(
                "Only a single group in the data. Can't run DE on a single group."
            )

    if not isinstance(group1, IterableClass) or isinstance(group1, str):
        group1 = [group1]

    # make a temp obs key using indices
    temp_key = None
    if idx1 is not None:
        obs_col, group1, group2 = _prepare_obs(idx1, idx2, adata)
        temp_key = "_scvi_temp_de"
        adata.obs[temp_key] = obs_col
        groupby = temp_key

    df_results = []
    dc = DifferentialComputation(model_fn, adata)
    for g1 in track(
        group1,
        description="DE...",
        disable=silent,
    ):
        cell_idx1 = (adata.obs[groupby] == g1).to_numpy().ravel()
        if group2 is None:
            cell_idx2 = ~cell_idx1
        else:
            cell_idx2 = (adata.obs[groupby] == group2).to_numpy().ravel()

        all_info = dc.get_bayes_factors(
            cell_idx1,
            cell_idx2,
            mode=mode,
            delta=delta,
            batchid1=batchid1,
            batchid2=batchid2,
            use_observed_batches=not batch_correction,
            **kwargs,
        )

        if all_stats is True:
            genes_properties_dict = all_stats_fn(adata, cell_idx1, cell_idx2)
            all_info = {**all_info, **genes_properties_dict}

        res = pd.DataFrame(all_info, index=col_names)
        sort_key = "proba_de" if mode == "change" else "bayes_factor"
        res = res.sort_values(by=sort_key, ascending=False)
        if mode == "change":
            res["is_de_fdr_{}".format(fdr)] = _fdr_de_prediction(
                res["proba_de"], fdr=fdr
            )
        if idx1 is None:
            g2 = "Rest" if group2 is None else group2
            res["comparison"] = "{} vs {}".format(g1, g2)
            res["group1"] = g1
            res["group2"] = g2
        df_results.append(res)

    if temp_key is not None:
        del adata.obs[temp_key]

    result = pd.concat(df_results, axis=0)

    return result


def _fdr_de_prediction(posterior_probas: np.ndarray, fdr: float = 0.05):
    """Compute posterior expected FDR and tag features as DE."""
    if not posterior_probas.ndim == 1:
        raise ValueError("posterior_probas should be 1-dimensional")
    sorted_genes = np.argsort(-posterior_probas)
    sorted_pgs = posterior_probas[sorted_genes]
    cumulative_fdr = (1.0 - sorted_pgs).cumsum() / (1.0 + np.arange(len(sorted_pgs)))
    d = (cumulative_fdr <= fdr).sum()
    pred_de_genes = sorted_genes[:d]
    is_pred_de = np.zeros_like(cumulative_fdr).astype(bool)
    is_pred_de[pred_de_genes] = True
    return is_pred_de
