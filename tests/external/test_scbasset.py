import numpy as np

from scvi.data import synthetic_iid
from scvi.external import SCBASSET

_DNA_CODE_KEY = "code"


def _get_adata(sparse=False):
    dataset1 = synthetic_iid(batch_size=100, sparse=sparse).transpose()
    dataset1.X = (dataset1.X > 0).astype(float)
    dataset1.obsm[_DNA_CODE_KEY] = np.random.randint(0, 3, size=(dataset1.n_obs, 1334))
    return dataset1


def test_scbasset():
    adata = _get_adata()
    SCBASSET.setup_anndata(
        adata,
        dna_code_key=_DNA_CODE_KEY,
    )
    model = SCBASSET(adata)
    model.train(max_epochs=2, early_stopping=True)
    model.get_latent_representation()


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
