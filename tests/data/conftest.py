import numpy as np
import pytest

from scvi.data import synthetic_iid


@pytest.fixture
def adata():
    anndata = synthetic_iid()
    raw_counts = anndata.X.copy()
    anndata.layers["raw"] = raw_counts
    anndata.obs["cont1"] = np.random.normal(size=(anndata.shape[0],))
    anndata.obs["cont2"] = np.random.normal(size=(anndata.shape[0],))
    anndata.obs["cat1"] = np.random.randint(0, 5, size=(anndata.shape[0],))
    anndata.obs["cat2"] = np.random.randint(0, 5, size=(anndata.shape[0],))

    return anndata


adata1 = adata2 = adata
