
from scvi.external import POISSONVI


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
