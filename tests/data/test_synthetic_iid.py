import numpy as np
import pytest

from scvi.data import synthetic_iid


@pytest.mark.parametrize("sparse_format", ["csr_matrix", "csc_matrix", None])
def test_synthetic_iid_sparse_format(sparse_format: str | None):
    _ = synthetic_iid(sparse_format=sparse_format)


@pytest.mark.parametrize("return_mudata", [False, True])
def test_synthetic_iid_coords(return_mudata: bool):
    adata = synthetic_iid(return_mudata=return_mudata, generate_coordinates=True)
    assert "coordinates" in adata.obsm
    assert isinstance(adata.obsm["coordinates"], np.ndarray)
    assert adata.obsm["coordinates"].shape == (adata.n_obs, 2)
