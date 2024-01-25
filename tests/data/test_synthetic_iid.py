import pytest

from scvi.data import synthetic_iid


@pytest.mark.parametrize("sparse_format", ["csr_matrix", "csc_matrix", None])
def test_synthetic_iid_sparse_format(sparse_format: str | None):
    _ = synthetic_iid(sparse_format=sparse_format)
