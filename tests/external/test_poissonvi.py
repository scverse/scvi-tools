from typing import Optional

from scvi.data import synthetic_iid
from scvi.external import POISSONVI


def _get_adata(sparse_format: Optional[str] = None):
    dataset1 = synthetic_iid(batch_size=100, sparse_format=sparse_format)
    return dataset1


def test_poissonvi():
    adata = _get_adata()
    print(adata.obs)
    POISSONVI.setup_anndata(
        adata,
    )
    model = POISSONVI(adata)
    model.train(max_epochs=1)
    model.get_latent_representation()
    model.get_accessibility_estimates()
