import logging
from typing import Optional, Union

import numpy as np

from scvi.data import AnnDataManager

logger = logging.getLogger(__name__)


def scmc_raw_counts_properties(
    adata_manager: AnnDataManager,
    idx1: Union[list[int], np.ndarray],
    idx2: Union[list[int], np.ndarray],
    var_idx: Optional[Union[list[int], np.ndarray]] = None,
    modality: str = None,
) -> dict[str, np.ndarray]:
    """Computes and returns some statistics on the raw counts of two sub-populations.

    Parameters
    ----------
    adata_manager
        :class:`~scvi.data.AnnDataManager` object setup with :class:`~scvi.model.SCVI`.
    idx1
        subset of indices describing the first population.
    idx2
        subset of indices describing the second population.
    mc_key
        denotes layer where methylated cytosine counts are stored in AnnData
    cov_key
        denotes layer where total cytosine counts are stored in AnnData
    var_idx
        subset of variables to extract properties from. if None, all variables are used.

    Returns
    -------
    type
        Dict of ``np.ndarray`` containing, by pair (one for each sub-population).
    """
    adata = adata_manager.adata
    mc = adata[modality].layers[adata_manager.mc_layer]
    cov = adata[modality].layers[adata_manager.cov_layer]
    mc1 = mc[idx1]
    mc2 = mc[idx2]

    cov1 = cov[idx1]
    cov2 = cov[idx2]
    if var_idx is not None:
        mc1 = mc1[:, var_idx]
        mc2 = mc2[:, var_idx]

        cov1 = cov1[:, var_idx]
        cov2 = cov2[:, var_idx]

    mean1 = np.asarray(np.nanmean(mc1 / cov1, axis=0)).ravel()
    mean2 = np.asarray(np.nanmean(mc2 / cov2, axis=0)).ravel()

    properties = {"emp_mean1": mean1, "emp_mean2": mean2, "emp_effect": (mean1 - mean2)}
    return properties
