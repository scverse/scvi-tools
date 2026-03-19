"""Tests for shared memory utilities."""

import os
import subprocess

import numpy as np
import pytest
import torch
from scipy.sparse import csr_matrix
from scipy.sparse import random as sparse_random

from scvi.dataloaders._shared_memory import (
    SharedMemoryRegistry,
    attach_dense_array,
    attach_sparse_csr,
    share_dense_array,
    share_sparse_csr,
)


@pytest.fixture
def dense_array():
    rng = np.random.default_rng(42)
    return rng.standard_normal((100, 50)).astype(np.float32)


@pytest.fixture
def sparse_csr():
    mat = sparse_random(100, 50, density=0.1, format="csr", dtype=np.float32, random_state=42)
    return mat


class TestShareDenseArray:
    def test_round_trip(self, dense_array):
        """Test creating and attaching shared memory for a dense array."""
        name = "test_dense_rt"
        shm, shared_arr = share_dense_array(name, dense_array)
        try:
            # Verify the shared array matches the original
            np.testing.assert_array_equal(shared_arr, dense_array)

            # Attach from a "different rank"
            shm2, attached_arr = attach_dense_array(name, dense_array.shape, dense_array.dtype)
            np.testing.assert_array_equal(attached_arr, dense_array)

            # Verify they share the same memory (modifying one affects the other)
            attached_arr[0, 0] = 999.0
            assert shared_arr[0, 0] == 999.0

            shm2.close()
        finally:
            shm.close()
            shm.unlink()

    def test_different_dtypes(self):
        """Test with various numpy dtypes."""
        for dtype in [np.float32, np.float64, np.int32, np.int64]:
            name = f"test_dtype_{dtype.__name__}"
            arr = np.arange(20, dtype=dtype).reshape(4, 5)
            shm, shared = share_dense_array(name, arr)
            try:
                shm2, attached = attach_dense_array(name, arr.shape, arr.dtype)
                np.testing.assert_array_equal(attached, arr)
                assert attached.dtype == dtype
                shm2.close()
            finally:
                shm.close()
                shm.unlink()


class TestShareSparseCsr:
    def test_round_trip(self, sparse_csr):
        """Test creating and attaching shared memory for a sparse CSR matrix."""
        name = "test_sparse_rt"
        metadata, shm_handles = share_sparse_csr(name, sparse_csr)
        try:
            # Attach from a "different rank"
            mat, shm_handles2 = attach_sparse_csr(name, metadata)

            # Verify the reconstructed matrix matches the original
            np.testing.assert_array_almost_equal(mat.toarray(), sparse_csr.toarray())
            assert mat.shape == sparse_csr.shape
            assert mat.nnz == sparse_csr.nnz

            for h in shm_handles2:
                h.close()
        finally:
            for h in shm_handles:
                h.close()
                h.unlink()

    def test_shared_data_modification(self, sparse_csr):
        """Test that modifying the attached sparse matrix affects the shared one."""
        name = "test_sparse_mod"
        metadata, shm_handles = share_sparse_csr(name, sparse_csr)
        try:
            mat, shm_handles2 = attach_sparse_csr(name, metadata)

            # Modify the data array via the attached matrix
            if mat.nnz > 0:
                mat.data[0] = 12345.0

                # Re-attach and verify the change is visible
                mat2, shm_handles3 = attach_sparse_csr(name, metadata)
                assert mat2.data[0] == 12345.0

                for h in shm_handles3:
                    h.close()

            for h in shm_handles2:
                h.close()
        finally:
            for h in shm_handles:
                h.close()
                h.unlink()

    def test_metadata_fields(self, sparse_csr):
        """Test that metadata contains the expected fields."""
        name = "test_sparse_meta"
        metadata, shm_handles = share_sparse_csr(name, sparse_csr)
        try:
            assert "shape" in metadata
            assert "data_dtype" in metadata
            assert "indices_dtype" in metadata
            assert "indptr_dtype" in metadata
            assert metadata["shape"] == sparse_csr.shape
        finally:
            for h in shm_handles:
                h.close()
                h.unlink()


class TestSharedMemoryRegistry:
    def test_register_and_cleanup(self):
        """Test registering handles and cleaning them up."""
        registry = SharedMemoryRegistry(is_rank0=True)

        # Create some shared memory blocks
        arr = np.zeros(100, dtype=np.float32)
        shm1, _ = share_dense_array("test_reg_1", arr)
        shm2, _ = share_dense_array("test_reg_2", arr)

        registry.register(shm1)
        registry.register(shm2)

        # Cleanup should close and unlink all
        registry.cleanup()

        # Verify they're unlinked (can't reattach)
        with pytest.raises(FileNotFoundError):
            attach_dense_array("test_reg_1", (100,), np.float32)
        with pytest.raises(FileNotFoundError):
            attach_dense_array("test_reg_2", (100,), np.float32)

    def test_register_list(self):
        """Test registering a list of handles at once."""
        registry = SharedMemoryRegistry(is_rank0=True)

        arr = np.zeros(100, dtype=np.float32)
        shm1, _ = share_dense_array("test_reglist_1", arr)
        shm2, _ = share_dense_array("test_reglist_2", arr)

        registry.register([shm1, shm2])
        registry.cleanup()

        with pytest.raises(FileNotFoundError):
            attach_dense_array("test_reglist_1", (100,), np.float32)

    def test_non_rank0_does_not_unlink(self):
        """Test that non-rank-0 registry closes but doesn't unlink."""
        # Rank 0 creates
        arr = np.zeros(100, dtype=np.float32)
        shm_creator, _ = share_dense_array("test_noreg_unlink", arr)

        # Non-rank-0 attaches
        shm_attacher, _ = attach_dense_array("test_noreg_unlink", (100,), np.float32)
        registry = SharedMemoryRegistry(is_rank0=False)
        registry.register(shm_attacher)
        registry.cleanup()

        # Should still be accessible (not unlinked)
        shm_check, _ = attach_dense_array("test_noreg_unlink", (100,), np.float32)
        shm_check.close()

        # Now clean up as rank 0
        shm_creator.close()
        shm_creator.unlink()

    def test_double_cleanup_is_safe(self):
        """Test that calling cleanup twice doesn't raise."""
        registry = SharedMemoryRegistry(is_rank0=True)
        arr = np.zeros(100, dtype=np.float32)
        shm, _ = share_dense_array("test_double_cleanup", arr)
        registry.register(shm)

        registry.cleanup()
        registry.cleanup()  # should be a no-op


class TestIsShareable:
    def test_dense_numpy(self):
        from scvi.dataloaders._shared_memory import _is_shareable

        assert _is_shareable(np.zeros((10, 5)))

    def test_sparse_scipy(self):
        from scvi.dataloaders._shared_memory import _is_shareable

        mat = csr_matrix(np.zeros((10, 5)))
        assert _is_shareable(mat)

    def test_h5py_dataset(self):
        import tempfile

        import h5py

        from scvi.dataloaders._shared_memory import _is_shareable

        with tempfile.NamedTemporaryFile(suffix=".h5") as f:
            with h5py.File(f.name, "w") as hf:
                hf.create_dataset("X", data=np.zeros((10, 5)))
            with h5py.File(f.name, "r") as hf:
                assert not _is_shareable(hf["X"])


def _launch_ddp_script(script_code: str, save_path: str, script_name: str):
    """Helper to write a script and launch it with torchrun."""
    temp_file_path = os.path.join(save_path, script_name)
    with open(temp_file_path, "w") as f:
        f.write(script_code)

    world_size = torch.cuda.device_count()
    command = [
        "torchrun",
        f"--nproc_per_node={world_size}",
        temp_file_path,
    ]
    try:
        subprocess.run(command, check=True)
    finally:
        if os.path.exists(temp_file_path):
            os.remove(temp_file_path)


@pytest.mark.multigpu
def test_shared_memory_ddp_dense(save_path: str):
    """Test shared memory with dense adata.X in DDP training."""
    training_code = """\
import torch
import scvi
from scvi.model import SCVI

adata = scvi.data.synthetic_iid()
SCVI.setup_anndata(adata)

model = SCVI(adata)
model.train(
    max_epochs=1,
    check_val_every_n_epoch=1,
    accelerator="gpu",
    devices=-1,
    strategy="ddp_find_unused_parameters_true",
    datasplitter_kwargs={"share_memory": True},
)
assert model.is_trained
"""
    _launch_ddp_script(training_code, save_path, "test_shm_dense_ddp.py")


@pytest.mark.multigpu
def test_shared_memory_ddp_sparse(save_path: str):
    """Test shared memory with sparse CSR adata.X in DDP training."""
    training_code = """\
import torch
import scipy.sparse
import scvi
from scvi.model import SCVI

adata = scvi.data.synthetic_iid()
# Convert X to sparse CSR
adata.X = scipy.sparse.csr_matrix(adata.X)
SCVI.setup_anndata(adata)

model = SCVI(adata)
model.train(
    max_epochs=1,
    check_val_every_n_epoch=1,
    accelerator="gpu",
    devices=-1,
    strategy="ddp_find_unused_parameters_true",
    datasplitter_kwargs={"share_memory": True},
)
assert model.is_trained
"""
    _launch_ddp_script(training_code, save_path, "test_shm_sparse_ddp.py")


@pytest.mark.multigpu
def test_shared_memory_ddp_disabled(save_path: str):
    """Test that training works with share_memory=False in DDP."""
    training_code = """\
import torch
import scvi
from scvi.model import SCVI

adata = scvi.data.synthetic_iid()
SCVI.setup_anndata(adata)

model = SCVI(adata)
model.train(
    max_epochs=1,
    check_val_every_n_epoch=1,
    accelerator="gpu",
    devices=-1,
    strategy="ddp_find_unused_parameters_true",
    datasplitter_kwargs={"share_memory": False},
)
assert model.is_trained
"""
    _launch_ddp_script(training_code, save_path, "test_shm_disabled_ddp.py")
