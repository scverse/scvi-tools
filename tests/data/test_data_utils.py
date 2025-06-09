import numpy as np
import pytest
import requests
import torch

import scvi
from scvi.data._utils import scipy_to_torch_sparse


@pytest.mark.parametrize("sparse_format", ["csr_matrix", "csc_matrix", None])
def test_scipy_to_torch_sparse(sparse_format: str):
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
            _ = scipy_to_torch_sparse(scipy_sparse)
        return

    torch_sparse = scipy_to_torch_sparse(scipy_sparse)
    assert isinstance(torch_sparse, torch.Tensor)
    assert torch_sparse.layout is expected_torch_layout
    assert torch_sparse.shape == scipy_sparse.shape
    assert np.allclose(
        torch_sparse.to_dense().numpy(),
        scipy_sparse.toarray(),
    )


def test_colab_link_alive():
    url = (
        "https://colab.research.google.com/github/scverse/scvi-tutorials/blob/main/"
        "quick_start/api_overview.ipynb"
    )
    response = requests.get(url, allow_redirects=True)
    assert response.status_code == 200, f"Link is not reachable, got {response.status_code}"
