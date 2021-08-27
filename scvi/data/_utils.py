import logging
from typing import Union

import anndata
import h5py
import numpy as np
import pandas as pd
import scipy.sparse as sp_sparse
from anndata._core.sparse_dataset import SparseDataset
from numba import boolean, float32, float64, int32, int64, vectorize

logger = logging.getLogger(__name__)


def _check_nonnegative_integers(
    data: Union[pd.DataFrame, np.ndarray, sp_sparse.spmatrix, h5py.Dataset]
):
    """Approximately checks values of data to ensure it is count data."""

    # for backed anndata
    if isinstance(data, h5py.Dataset) or isinstance(data, SparseDataset):
        data = data[:100]

    if isinstance(data, np.ndarray):
        data = data
    elif issubclass(type(data), sp_sparse.spmatrix):
        data = data.data
    elif isinstance(data, pd.DataFrame):
        data = data.to_numpy()
    else:
        raise TypeError("data type not understood")

    n = len(data)
    inds = np.random.permutation(n)[:20]
    check = data.flat[inds]
    return ~np.any(_is_not_count(check))


@vectorize(
    [
        boolean(int32),
        boolean(int64),
        boolean(float32),
        boolean(float64),
    ],
    target="parallel",
    cache=True,
)
def _is_not_count(d):
    return d < 0 or d % 1 != 0


def _get_batch_mask_protein_data(
    adata: anndata.AnnData, protein_expression_obsm_key: str, batch_key: str
):
    """
    Returns a list with length number of batches where each entry is a mask.

    The mask is over cell measurement columns that are present (observed)
    in each batch. Absence is defined by all 0 for that protein in that batch.
    """
    pro_exp = adata.obsm[protein_expression_obsm_key]
    pro_exp = pro_exp.to_numpy() if isinstance(pro_exp, pd.DataFrame) else pro_exp
    batches = adata.obs[batch_key].values
    batch_mask = {}
    for b in np.unique(batches):
        b_inds = np.where(batches.ravel() == b)[0]
        batch_sum = pro_exp[b_inds, :].sum(axis=0)
        all_zero = batch_sum == 0
        batch_mask[b] = ~all_zero

    return batch_mask
