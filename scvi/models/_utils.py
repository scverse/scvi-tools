import anndata
import numpy as np
from scvi.dataset import get_from_registry
from typing import Union, Tuple, List
from scvi import _CONSTANTS
import logging

logger = logging.getLogger(__name__)


def scrna_raw_counts_properties(
    adata: anndata.AnnData,
    idx1: Union[List[int], np.ndarray],
    idx2: Union[List[int], np.ndarray],
) -> Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """
    Computes and returns some statistics on the raw counts of two sub-populations.

    Parameters
    ----------
    adata
        AnnData object setup with `scvi`.
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
    data = get_from_registry(adata, _CONSTANTS.X_KEY)
    data1 = data[idx1]
    data2 = data[idx2]
    mean1 = np.asarray((data1).mean(axis=0)).ravel()
    mean2 = np.asarray((data2).mean(axis=0)).ravel()
    nonz1 = np.asarray((data1 != 0).mean(axis=0)).ravel()
    nonz2 = np.asarray((data2 != 0).mean(axis=0)).ravel()

    key = "_scvi_raw_norm_X"
    if key not in adata.obsm.keys():
        scaling_factor = np.asarray(data.sum(axis=1)).ravel().reshape(-1, 1)
        normalized_data = 1e4 * data / scaling_factor
        adata.obsm[key] = normalized_data
        logger.info(
            "Storing library size normalized data in adata.obsm['{}']".format(key)
        )
        logger.info(
            "This can deleted after DE inference with `del adata.obsm['{}']`".format(
                key
            )
        )
    else:
        normalized_data = adata.obsm[key]

    norm_mean1 = np.asarray(normalized_data[idx1, :].mean(axis=0)).ravel()
    norm_mean2 = np.asarray(normalized_data[idx2, :].mean(axis=0)).ravel()
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


def _get_batch_code_from_category(adata, category):
    categorical_mappings = adata.uns["_scvi"]["categorical_mappings"]
    batch_mappings = categorical_mappings["_scvi_batch"]["mapping"]
    if category not in batch_mappings:
        raise ValueError('"{}" not a valid batch category.'.format(category))
    return np.where(batch_mappings == category)[0][0]
