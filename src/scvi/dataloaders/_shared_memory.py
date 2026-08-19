"""Shared memory utilities for DDP data deduplication.

In subprocess-based DDP, each GPU rank re-executes the full script and loads
its own copy of ``adata`` into memory. This module provides utilities to share
a single physical copy of ``adata.X`` across all ranks on the same node using
POSIX shared memory.
"""

from __future__ import annotations

import atexit
import gc
import logging
import os
import platform
from multiprocessing import shared_memory

import numpy as np
from scipy.sparse import csr_matrix, issparse

logger = logging.getLogger(__name__)


def share_dense_array(
    name: str,
    arr: np.ndarray,
) -> tuple[shared_memory.SharedMemory, np.ndarray]:
    """Create shared memory from a dense numpy array (rank 0).

    Parameters
    ----------
    name
        Name for the shared memory block.
    arr
        Dense numpy array to share.

    Returns
    -------
    Tuple of (SharedMemory handle, numpy view into shared memory).
    """
    shm = shared_memory.SharedMemory(name=name, create=True, size=arr.nbytes)
    shared_arr = np.ndarray(arr.shape, dtype=arr.dtype, buffer=shm.buf)
    np.copyto(shared_arr, arr)
    return shm, shared_arr


def attach_dense_array(
    name: str,
    shape: tuple[int, ...],
    dtype: np.dtype,
) -> tuple[shared_memory.SharedMemory, np.ndarray]:
    """Attach to existing shared memory as a dense numpy array (rank 1+).

    Parameters
    ----------
    name
        Name of the existing shared memory block.
    shape
        Shape of the array.
    dtype
        Data type of the array.

    Returns
    -------
    Tuple of (SharedMemory handle, numpy view into shared memory).
    """
    shm = shared_memory.SharedMemory(name=name, create=False)
    shared_arr = np.ndarray(shape, dtype=dtype, buffer=shm.buf)
    return shm, shared_arr


def share_sparse_csr(
    name: str,
    mat: csr_matrix,
) -> tuple[dict, list[shared_memory.SharedMemory]]:
    """Create shared memory from a sparse CSR matrix (rank 0).

    Shares 3 arrays (data, indices, indptr) as separate shared memory blocks.

    Parameters
    ----------
    name
        Base name for the shared memory blocks.
    mat
        Sparse CSR matrix to share.

    Returns
    -------
    Tuple of (metadata dict for broadcast, list of SharedMemory handles).
    """
    mat = csr_matrix(mat)  # ensure CSR format

    shm_data, shared_data = share_dense_array(f"{name}_data", mat.data)
    shm_indices, shared_indices = share_dense_array(f"{name}_indices", mat.indices)
    shm_indptr, shared_indptr = share_dense_array(f"{name}_indptr", mat.indptr)

    metadata = {
        "shape": mat.shape,
        "data_dtype": str(mat.data.dtype),
        "data_shape": mat.data.shape,
        "indices_dtype": str(mat.indices.dtype),
        "indices_shape": mat.indices.shape,
        "indptr_dtype": str(mat.indptr.dtype),
        "indptr_shape": mat.indptr.shape,
    }
    return metadata, [shm_data, shm_indices, shm_indptr]


def attach_sparse_csr(
    name: str,
    metadata: dict,
) -> tuple[csr_matrix, list[shared_memory.SharedMemory]]:
    """Attach to existing shared memory and reconstruct a sparse CSR matrix (rank 1+).

    Parameters
    ----------
    name
        Base name for the shared memory blocks.
    metadata
        Metadata dict from :func:`share_sparse_csr` (broadcast from rank 0).

    Returns
    -------
    Tuple of (CSR matrix backed by shared memory, list of SharedMemory handles).
    """
    shm_data, shared_data = attach_dense_array(
        f"{name}_data",
        metadata["data_shape"],
        np.dtype(metadata["data_dtype"]),
    )
    shm_indices, shared_indices = attach_dense_array(
        f"{name}_indices",
        metadata["indices_shape"],
        np.dtype(metadata["indices_dtype"]),
    )
    shm_indptr, shared_indptr = attach_dense_array(
        f"{name}_indptr",
        metadata["indptr_shape"],
        np.dtype(metadata["indptr_dtype"]),
    )

    mat = csr_matrix(
        (shared_data, shared_indices, shared_indptr),
        shape=metadata["shape"],
        copy=False,
    )
    return mat, [shm_data, shm_indices, shm_indptr]


class SharedMemoryRegistry:
    """Tracks shared memory blocks and handles cleanup.

    Parameters
    ----------
    is_rank0
        Whether this process is rank 0 (responsible for unlinking).
    """

    def __init__(self, is_rank0: bool = False):
        self.is_rank0 = is_rank0
        self._handles: list[shared_memory.SharedMemory] = []
        self._cleaned_up = False

    def register(self, shm: shared_memory.SharedMemory | list[shared_memory.SharedMemory]):
        """Register shared memory handle(s) for cleanup."""
        if isinstance(shm, list):
            self._handles.extend(shm)
        else:
            self._handles.append(shm)

    def cleanup(self):
        """Close all handles and unlink on rank 0."""
        if self._cleaned_up:
            return
        self._cleaned_up = True

        for shm in self._handles:
            try:
                shm.close()
            except OSError:
                pass
            if self.is_rank0:
                try:
                    shm.unlink()
                except OSError:
                    pass
        self._handles.clear()


def _malloc_trim():
    """Call malloc_trim on Linux to return freed pages to the OS."""
    if platform.system() == "Linux":
        try:
            import ctypes

            libc = ctypes.CDLL("libc.so.6")
            libc.malloc_trim(0)
        except OSError:
            pass


def _is_shareable(X) -> bool:
    """Check if the data matrix X is shareable (dense numpy or scipy sparse).

    Returns False for backed (h5py) or dask arrays.
    """
    import h5py

    if isinstance(X, h5py.Dataset):
        return False
    if isinstance(X, np.ndarray):
        return True
    if issparse(X):
        return True
    try:
        import dask.array as da

        if isinstance(da.Array, type) and isinstance(X, da.Array):
            return False
    except ImportError:
        pass
    return False


def setup_shared_memory(
    adata_manager,
    registry: SharedMemoryRegistry | None = None,
) -> SharedMemoryRegistry | None:
    """Orchestrate shared memory setup for DDP.

    Called from ``DataSplitter.setup()`` after index splitting. Rank 0 copies
    ``adata.X`` into shared memory; all ranks replace ``adata.X`` with a view
    into the shared block.

    Parameters
    ----------
    adata_manager
        The AnnDataManager whose adata.X should be shared.
    registry
        Existing registry, or None to create a new one.

    Returns
    -------
    SharedMemoryRegistry if shared memory was set up, None otherwise.
    """
    import torch.distributed as dist

    if not dist.is_initialized():
        return None

    rank = dist.get_rank()
    world_size = dist.get_world_size()
    if world_size <= 1:
        return None

    adata = adata_manager.adata
    X = adata.X

    if not _is_shareable(X):
        logger.info("adata.X is not shareable (backed or dask); skipping shared memory.")
        return None

    if registry is None:
        registry = SharedMemoryRegistry(is_rank0=(rank == 0))

    # Use rank 0's PID for unique naming across training sessions
    pid_list = [os.getpid() if rank == 0 else 0]
    dist.broadcast_object_list(pid_list, src=0)
    base_name = f"scvi_{pid_list[0]}_X"

    is_sparse = issparse(X)

    if rank == 0:
        # Clean up stale shared memory with the same name
        _cleanup_stale_shm(base_name, is_sparse)

        if is_sparse:
            metadata, shm_handles = share_sparse_csr(base_name, X)
            registry.register(shm_handles)
            broadcast_data = [{"sparse": True, "metadata": metadata}]
        else:
            shm, shared_arr = share_dense_array(base_name, X)
            registry.register(shm)
            broadcast_data = [
                {
                    "sparse": False,
                    "shape": X.shape,
                    "dtype": str(X.dtype),
                }
            ]
        logger.info(
            f"[rank 0] Created shared memory for adata.X "
            f"({X.nbytes / 1024**2:.1f} MB, sparse={is_sparse})"
        )
    else:
        broadcast_data = [None]

    dist.broadcast_object_list(broadcast_data, src=0)
    info = broadcast_data[0]

    if rank != 0:
        if info["sparse"]:
            shared_X, shm_handles = attach_sparse_csr(base_name, info["metadata"])
            registry.register(shm_handles)
        else:
            shm, shared_X = attach_dense_array(base_name, info["shape"], np.dtype(info["dtype"]))
            registry.register(shm)
        logger.info(f"[rank {rank}] Attached to shared memory for adata.X")
    else:
        if is_sparse:
            # Rank 0: reconstruct from shared memory too so we can free the original
            shared_X, attach_handles = attach_sparse_csr(base_name, info["metadata"])
            registry.register(attach_handles)
        else:
            shared_X = shared_arr

    # Replace adata.X with the shared view
    del X
    adata.X = shared_X
    gc.collect()
    _malloc_trim()

    # Register atexit cleanup as a safety net
    atexit.register(registry.cleanup)

    dist.barrier()
    return registry


def _cleanup_stale_shm(base_name: str, is_sparse: bool):
    """Try to clean up stale shared memory from a previous run."""
    names = (
        [base_name]
        if not is_sparse
        else [f"{base_name}_data", f"{base_name}_indices", f"{base_name}_indptr"]
    )
    for name in names:
        try:
            old = shared_memory.SharedMemory(name=name, create=False)
            old.close()
            old.unlink()
        except FileNotFoundError:
            pass
