import anndata
import numpy as np
from scvi.dataset import get_from_registry
from typing import Union, Tuple, List
from scvi import _CONSTANTS


def scrna_raw_counts_properties(
    adata: anndata.AnnData,
    idx1: Union[List[int], np.ndarray],
    idx2: Union[List[int], np.ndarray],
) -> Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """
    Computes and returns some statistics on the raw counts of two sub-populations.

    Parameters
    ----------
    idx1
        subset of indices describing the first population.
    idx2
        subset of indices describing the second population.
    Returns
    -------
    type
        Tuple of ``np.ndarray`` containing, by pair (one for each sub-population),
        mean expression per gene, proportion of non-zero expression per gene, mean of normalized expression.
    """
    X = get_from_registry(adata, _CONSTANTS.X_KEY)
    mean1 = np.asarray((X[idx1]).mean(axis=0)).ravel()
    mean2 = np.asarray((X[idx2]).mean(axis=0)).ravel()
    nonz1 = np.asarray((X[idx1] != 0).mean(axis=0)).ravel()
    nonz2 = np.asarray((X[idx2] != 0).mean(axis=0)).ravel()

    scaling_factor = np.asarray(X.sum(axis=1)).ravel().reshape(-1, 1)
    normalized_X = X / scaling_factor

    norm_mean1 = np.asarray(normalized_X[idx1, :].mean(axis=0)).ravel()
    norm_mean2 = np.asarray(normalized_X[idx2, :].mean(axis=0)).ravel()
    return_vals = [mean1, mean2, nonz1, nonz2, norm_mean1, norm_mean2]
    return [np.squeeze(np.asarray(arr)) for arr in return_vals]


def _get_var_names_from_setup_anndata(adata):
    """Gets var names by checking if using raw."""
    var_names = (
        adata.var_names
        if adata.uns["_scvi"]["use_raw"] is False
        else adata.raw.var_names
    )

    return var_names
