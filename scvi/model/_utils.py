import logging
import warnings
from collections.abc import Iterable as IterableClass
from typing import Dict, List, Optional, Sequence, Tuple, Union

import numpy as np
import scipy.sparse as sp_sparse
import torch

from scvi import REGISTRY_KEYS
from scvi._types import Number
from scvi.data import AnnDataManager

logger = logging.getLogger(__name__)


def parse_use_gpu_arg(
    use_gpu: Optional[Union[str, int, bool]] = None,
    return_device=True,
):
    """
    Parses the use_gpu arg in codebase.

    Returned gpus are is compatible with PytorchLightning's gpus arg.
    If return_device is True, will also return the device.

    Parameters
    ----------
    use_gpu
        Use default GPU if available (if None or True), or index of GPU to use (if int),
        or name of GPU (if str, e.g., `'cuda:0'`), or use CPU (if False).
    return_device
        If True, will return the torch.device of use_gpu.

    Returns
    -------
    Arguments for lightning trainer, including the accelerator (str), devices
    (int or sequence of int), and optionally the torch device.
    """
    # Support Apple silicon
    cuda_available = torch.cuda.is_available()
    # If using an older version of torch.
    try:
        mps_available = torch.backends.mps.is_available()
    except AttributeError:
        mps_available = False
    gpu_available = cuda_available
    lightning_devices = None
    if (use_gpu is None and not gpu_available) or (use_gpu is False):
        accelerator = "cpu"
        device = torch.device("cpu")
    elif (use_gpu is None and gpu_available) or (use_gpu is True):
        current = torch.cuda.current_device() if cuda_available else "mps"
        if current != "mps":
            lightning_devices = [current]
            accelerator = "gpu"
        else:
            accelerator = "mps"
            lightning_devices = 1
        device = torch.device(current)
    # Also captures bool case
    elif isinstance(use_gpu, int):
        device = torch.device(use_gpu) if not mps_available else torch.device("mps")
        accelerator = "gpu" if not mps_available else "mps"
        lightning_devices = [use_gpu] if not mps_available else 1
    elif isinstance(use_gpu, str):
        device = torch.device(use_gpu)
        accelerator = "gpu"
        # changes "cuda:0" to "0,"
        lightning_devices = [int(use_gpu.split(":")[-1])]
    else:
        raise ValueError("use_gpu argument not understood.")

    if return_device:
        return accelerator, lightning_devices, device
    else:
        return accelerator, lightning_devices


def scrna_raw_counts_properties(
    adata_manager: AnnDataManager,
    idx1: Union[List[int], np.ndarray],
    idx2: Union[List[int], np.ndarray],
    var_idx: Optional[Union[List[int], np.ndarray]] = None,
) -> Dict[str, np.ndarray]:
    """
    Computes and returns some statistics on the raw counts of two sub-populations.

    Parameters
    ----------
    adata_manager
        :class:`~scvi.data.AnnDataManager` object setup with :class:`~scvi.model.SCVI`.
    idx1
        subset of indices describing the first population.
    idx2
        subset of indices describing the second population.
    var_idx
        subset of variables to extract properties from. if None, all variables are used.

    Returns
    -------
    type
        Dict of ``np.ndarray`` containing, by pair (one for each sub-population),
        mean expression per gene, proportion of non-zero expression per gene, mean of normalized expression.
    """
    adata = adata_manager.adata
    data = adata_manager.get_from_registry(REGISTRY_KEYS.X_KEY)
    data1 = data[idx1]
    data2 = data[idx2]
    if var_idx is not None:
        data1 = data1[:, var_idx]
        data2 = data2[:, var_idx]

    mean1 = np.asarray((data1).mean(axis=0)).ravel()
    mean2 = np.asarray((data2).mean(axis=0)).ravel()
    nonz1 = np.asarray((data1 != 0).mean(axis=0)).ravel()
    nonz2 = np.asarray((data2 != 0).mean(axis=0)).ravel()

    key = "_scvi_raw_norm_scaling"
    if key not in adata.obs.keys():
        scaling_factor = 1 / np.asarray(data.sum(axis=1)).ravel().reshape(-1, 1)
        scaling_factor *= 1e4
        adata.obs[key] = scaling_factor.ravel()
    else:
        scaling_factor = adata.obs[key].to_numpy().ravel().reshape(-1, 1)

    if issubclass(type(data), sp_sparse.spmatrix):
        norm_data1 = data1.multiply(scaling_factor[idx1])
        norm_data2 = data2.multiply(scaling_factor[idx2])
    else:
        norm_data1 = data1 * scaling_factor[idx1]
        norm_data2 = data2 * scaling_factor[idx2]

    norm_mean1 = np.asarray(norm_data1.mean(axis=0)).ravel()
    norm_mean2 = np.asarray(norm_data2.mean(axis=0)).ravel()

    properties = dict(
        raw_mean1=mean1,
        raw_mean2=mean2,
        non_zeros_proportion1=nonz1,
        non_zeros_proportion2=nonz2,
        raw_normalized_mean1=norm_mean1,
        raw_normalized_mean2=norm_mean2,
    )
    return properties


def cite_seq_raw_counts_properties(
    adata_manager: AnnDataManager,
    idx1: Union[List[int], np.ndarray],
    idx2: Union[List[int], np.ndarray],
) -> Dict[str, np.ndarray]:
    """
    Computes and returns some statistics on the raw counts of two sub-populations.

    Parameters
    ----------
    adata_manager
        :class:`~scvi.data.AnnDataManager` object setup with :class:`~scvi.model.TOTALVI`.
    idx1
        subset of indices describing the first population.
    idx2
        subset of indices describing the second population.

    Returns
    -------
    type
        Dict of ``np.ndarray`` containing, by pair (one for each sub-population),
        mean expression per gene, proportion of non-zero expression per gene, mean of normalized expression.
    """
    gp = scrna_raw_counts_properties(adata_manager, idx1, idx2)
    protein_exp = adata_manager.get_from_registry(REGISTRY_KEYS.PROTEIN_EXP_KEY)

    nan = np.array([np.nan] * adata_manager.summary_stats.n_proteins)
    protein_exp = adata_manager.get_from_registry(REGISTRY_KEYS.PROTEIN_EXP_KEY)
    mean1_pro = np.asarray(protein_exp[idx1].mean(0))
    mean2_pro = np.asarray(protein_exp[idx2].mean(0))
    nonz1_pro = np.asarray((protein_exp[idx1] > 0).mean(0))
    nonz2_pro = np.asarray((protein_exp[idx2] > 0).mean(0))
    properties = dict(
        raw_mean1=np.concatenate([gp["raw_mean1"], mean1_pro]),
        raw_mean2=np.concatenate([gp["raw_mean2"], mean2_pro]),
        non_zeros_proportion1=np.concatenate([gp["non_zeros_proportion1"], nonz1_pro]),
        non_zeros_proportion2=np.concatenate([gp["non_zeros_proportion2"], nonz2_pro]),
        raw_normalized_mean1=np.concatenate([gp["raw_normalized_mean1"], nan]),
        raw_normalized_mean2=np.concatenate([gp["raw_normalized_mean2"], nan]),
    )

    return properties


def scatac_raw_counts_properties(
    adata_manager: AnnDataManager,
    idx1: Union[List[int], np.ndarray],
    idx2: Union[List[int], np.ndarray],
    var_idx: Optional[Union[List[int], np.ndarray]] = None,
) -> Dict[str, np.ndarray]:
    """
    Computes and returns some statistics on the raw counts of two sub-populations.

    Parameters
    ----------
    adata_manager
        :class:`~scvi.data.AnnDataManager` object setup with :class:`~scvi.model.SCVI`.
    idx1
        subset of indices describing the first population.
    idx2
        subset of indices describing the second population.
    var_idx
        subset of variables to extract properties from. if None, all variables are used.

    Returns
    -------
    type
        Dict of ``np.ndarray`` containing, by pair (one for each sub-population).
    """
    data = adata_manager.get_from_registry(REGISTRY_KEYS.X_KEY)
    data1 = data[idx1]
    data2 = data[idx2]
    if var_idx is not None:
        data1 = data1[:, var_idx]
        data2 = data2[:, var_idx]
    mean1 = np.asarray((data1 > 0).mean(axis=0)).ravel()
    mean2 = np.asarray((data2 > 0).mean(axis=0)).ravel()
    properties = dict(emp_mean1=mean1, emp_mean2=mean2, emp_effect=(mean1 - mean2))
    return properties


def _get_batch_code_from_category(
    adata_manager: AnnDataManager, category: Sequence[Union[Number, str]]
):
    if not isinstance(category, IterableClass) or isinstance(category, str):
        category = [category]

    batch_mappings = adata_manager.get_state_registry(
        REGISTRY_KEYS.BATCH_KEY
    ).categorical_mapping
    batch_code = []
    for cat in category:
        if cat is None:
            batch_code.append(None)
        elif cat not in batch_mappings:
            raise ValueError(f'"{cat}" not a valid batch category.')
        else:
            batch_loc = np.where(batch_mappings == cat)[0][0]
            batch_code.append(batch_loc)
    return batch_code


def _init_library_size(
    adata_manager: AnnDataManager, n_batch: dict
) -> Tuple[np.ndarray, np.ndarray]:
    """
    Computes and returns library size.

    Parameters
    ----------
    adata_manager
        :class:`~scvi.data.AnnDataManager` object setup with :class:`~scvi.model.SCVI`.
    n_batch
        Number of batches.

    Returns
    -------
    type
        Tuple of two 1 x n_batch ``np.ndarray`` containing the means and variances
        of library size in each batch in adata.

        If a certain batch is not present in the adata, the mean defaults to 0,
        and the variance defaults to 1. These defaults are arbitrary placeholders which
        should not be used in any downstream computation.
    """
    data = adata_manager.get_from_registry(REGISTRY_KEYS.X_KEY)
    batch_indices = adata_manager.get_from_registry(REGISTRY_KEYS.BATCH_KEY)

    library_log_means = np.zeros(n_batch)
    library_log_vars = np.ones(n_batch)

    for i_batch in np.unique(batch_indices):
        idx_batch = np.squeeze(batch_indices == i_batch)
        batch_data = data[
            idx_batch.nonzero()[0]
        ]  # h5ad requires integer indexing arrays.
        sum_counts = batch_data.sum(axis=1)
        masked_log_sum = np.ma.log(sum_counts)
        if np.ma.is_masked(masked_log_sum):
            warnings.warn(
                "This dataset has some empty cells, this might fail inference."
                "Data should be filtered with `scanpy.pp.filter_cells()`"
            )

        log_counts = masked_log_sum.filled(0)
        library_log_means[i_batch] = np.mean(log_counts).astype(np.float32)
        library_log_vars[i_batch] = np.var(log_counts).astype(np.float32)

    return library_log_means.reshape(1, -1), library_log_vars.reshape(1, -1)


def _get_var_names_from_manager(
    adata_manager: AnnDataManager, registry_key: str = REGISTRY_KEYS.X_KEY
) -> np.ndarray:
    return np.asarray(adata_manager.get_state_registry(registry_key).column_names)
