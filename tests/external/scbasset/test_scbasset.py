import numpy as np
import pytest

from scvi.data import synthetic_iid
from scvi.external import SCBASSET

_DNA_CODE_KEY = "code"


def _get_adata(sparse_format: str | None = None):
    dataset1 = synthetic_iid(batch_size=100, sparse_format=sparse_format).transpose()
    dataset1.X = (dataset1.X > 0).astype(float)
    dataset1.obsm[_DNA_CODE_KEY] = np.random.randint(0, 3, size=(dataset1.n_obs, 1344))
    return dataset1


@pytest.mark.optional
def test_scbasset():
    adata = _get_adata()
    SCBASSET.setup_anndata(
        adata,
        dna_code_key=_DNA_CODE_KEY,
    )
    model = SCBASSET(adata)
    model.train(max_epochs=2, early_stopping=True)
    model.get_latent_representation()


@pytest.mark.optional
def test_scbasset_batch():
    adata = _get_adata()
    SCBASSET.setup_anndata(
        adata,
        dna_code_key=_DNA_CODE_KEY,
        batch_key="batch",
    )
    model = SCBASSET(adata)
    model.train(max_epochs=1)
    model.get_latent_representation()
    assert hasattr(model.module, "batch_ids")
