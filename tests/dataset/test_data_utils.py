import numpy as np
import pytest
import torch

import scvi
from scvi.data._utils import convert_scipy_sparse_to_torch_sparse


@pytest.mark.parametrize(
    "sparse_format", ["csr_matrix", "csr_array", "csc_matrix", "csc_array", None]
)
def test_convert_scipy_sparse_to_torch_sparse(sparse_format: str):
    adata = scvi.data.synthetic_iid(sparse_format=sparse_format)
    scipy_sparse = adata.X

    if sparse_format is not None:
        if "csr" in sparse_format:
            expected_torch_layout = torch.sparse_csr
        else:
            expected_torch_layout = torch.sparse_csc
    else:
        # raises error for dense data
        with pytest.raises(TypeError):
            _ = convert_scipy_sparse_to_torch_sparse(scipy_sparse)
        return

    torch_sparse = convert_scipy_sparse_to_torch_sparse(scipy_sparse)
    assert isinstance(torch_sparse, torch.Tensor)
    assert torch_sparse.layout is expected_torch_layout
    assert torch_sparse.shape == scipy_sparse.shape
    assert np.allclose(
        torch_sparse.to_dense().numpy(),
        scipy_sparse.toarray(),
    )
