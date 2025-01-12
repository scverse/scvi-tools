import numpy as np
import pytest

from scvi.data import synthetic_iid
from scvi.external import RESOLVI


@pytest.fixture(scope="session")
def adata():
    adata = synthetic_iid(
        generate_coordinates=True,
        n_regions=0,
        n_proteins=0,
    )
    adata.obsm["X_spatial"] = adata.obsm["coordinates"]
    print(adata)
    return adata


def test_resolvi_train(adata):
    RESOLVI.setup_anndata(adata)
    model = RESOLVI(adata)
    model.train(
        max_epochs=2,
    )


def test_resolvi_downstream(adata):
    RESOLVI.setup_anndata(adata)
    model = RESOLVI(adata)
    model.train(
        max_epochs=2,
    )
    latent = model.get_latent_representation()
    assert latent.shape == (adata.n_obs, model.module.n_latent)
    model.differential_expression(groupby="labels")
    model.differential_expression(groupby="labels", weights="importance")
    model.save("test_resolvi")
    model2 = model.load("test_resolvi")
    latent2 = model2.get_latent_representation()
    assert np.allclose(latent, latent2)

    model_query = model.load_query_data(adata)
    model_query.train(
        max_epochs=2,
    )

def test_resolvi_semisupervised(adata):
    RESOLVI.setup_anndata(adata, label_key="labels")
    model = RESOLVI(adata, semisupervised=True)
    model.train(
        max_epochs=2,
        check_val_every_n_epoch=1,
        train_size=0.5,
        early_stopping=True,
    )
    model.differential_niche_abundance()
    model.predict()