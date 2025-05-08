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
    model = RESOLVI(adata, dispersion="gene-batch")
    model.train(
        max_epochs=2,
    )


def test_resolvi_save_load(adata):
    RESOLVI.setup_anndata(adata)
    model = RESOLVI(adata)
    model.train(
        max_epochs=2,
    )
    hist_elbo = model.history_["elbo_train"]
    latent = model.get_latent_representation()
    assert latent.shape == (adata.n_obs, model.module.n_latent)
    model.differential_expression(groupby="labels")
    model.differential_expression(groupby="labels", weights="importance")
    model.save("test_resolvi", save_anndata=True, overwrite=True)
    model2 = model.load("test_resolvi")
    np.testing.assert_array_equal(model2.history_["elbo_train"], hist_elbo)
    latent2 = model2.get_latent_representation()
    assert np.allclose(latent, latent2)
    model.load_query_data(reference_model="test_resolvi", adata=adata)


def test_resolvi_downstream(adata):
    RESOLVI.setup_anndata(adata)
    model = RESOLVI(adata)
    model.train(
        max_epochs=2,
    )
    latent = model.get_latent_representation()
    assert latent.shape == (adata.n_obs, model.module.n_latent)
    counts = model.get_normalized_expression(n_samples=31, library_size=10000)
    counts = model.get_normalized_expression_importance(n_samples=30, library_size=10000)
    print("FFFFFF", counts.shape)
    model.differential_expression(groupby="labels")
    model.differential_expression(groupby="labels", weights="importance")
    model.sample_posterior(
        model=model.module.model_residuals,
        num_samples=30,
        return_samples=False,
        return_sites=None,
        batch_size=1000,
    )
    model.sample_posterior(
        model=model.module.model_residuals, num_samples=30, return_samples=False, batch_size=1000
    )
    model_query = model.load_query_data(reference_model=model, adata=adata)
    model_query = model.load_query_data(reference_model="test_resolvi", adata=adata)
    model_query.train(
        max_epochs=2,
    )


def test_resolvi_semisupervised(adata):
    RESOLVI.setup_anndata(adata, labels_key="labels")
    model = RESOLVI(adata, semisupervised=True)
    model.train(
        max_epochs=2,
    )
    model.differential_niche_abundance(
        batch_size=30,
        groupby="batch",
        neighbor_key="index_neighbor",
    )
    pred = model.predict(soft=True)
    assert pred.shape == (adata.n_obs, model.summary_stats.n_labels - 1)
    pred = model.predict(soft=False)
    assert pred.shape == (adata.n_obs,)
