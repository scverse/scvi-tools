from scvi.data import synthetic_iid
from scvi.external import POISSONVI


def test_poissonvi():
    adata = synthetic_iid(batch_size=100)
    POISSONVI.setup_anndata(
        adata,
    )
    model = POISSONVI(adata)
    model.train(max_epochs=1)
    model.get_latent_representation()
    model.get_accessibility_estimates()
