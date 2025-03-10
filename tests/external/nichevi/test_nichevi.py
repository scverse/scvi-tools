import numpy as np
import pytest
from anndata import AnnData

from scvi.data import synthetic_iid
from scvi.external import nicheSCVI

N_LATENT = 10
K_NN = 5
N_EPOCHS_NICHEVI = 2

setup_kwargs = {
    "sample_key": "batch",
    "labels_key": "labels",
    "cell_coordinates_key": "coordinates",
    "expression_embedding_key": "qz1_m",
    "expression_embedding_niche_key": "qz1_m_niche_ct",
    "niche_composition_key": "neighborhood_composition",
    "niche_indexes_key": "niche_indexes",
    "niche_distances_key": "niche_distances",
}


@pytest.fixture(scope="session")
def adata():
    adata = synthetic_iid(
        batch_size=256,
        n_genes=100,
        n_proteins=0,
        n_regions=0,
        n_batches=2,
        n_labels=3,
        dropout_ratio=0.5,
        generate_coordinates=True,
        sparse_format=None,
        return_mudata=False,
    )

    adata.obsm["qz1_m"] = np.random.normal(size=(adata.shape[0], N_LATENT))
    adata.layers["counts"] = adata.X.copy()

    return adata


def test_nichevi_train(adata: AnnData):
    nicheSCVI.preprocessing_anndata(
        adata,
        k_nn=K_NN,
        **setup_kwargs,
    )

    nicheSCVI.setup_anndata(
        adata,
        layer="counts",
        batch_key="batch",
        **setup_kwargs,
    )
    nichevae = nicheSCVI(
        adata,
        prior_mixture=False,
        semisupervised=True,
        linear_classifier=True,
    )

    nichevae.train(
        max_epochs=N_EPOCHS_NICHEVI,
        train_size=0.8,
        validation_size=0.2,
        early_stopping=True,
        check_val_every_n_epoch=1,
        accelerator="cpu",
    )


def test_nichevi_save_load(adata):
    nicheSCVI.preprocessing_anndata(
        adata,
        k_nn=K_NN,
        **setup_kwargs,
    )

    nicheSCVI.setup_anndata(
        adata,
        layer="counts",
        batch_key="batch",
        **setup_kwargs,
    )
    nichevae = nicheSCVI(
        adata,
        prior_mixture=False,
        semisupervised=True,
        linear_classifier=True,
    )

    nichevae.train(
        max_epochs=N_EPOCHS_NICHEVI,
        train_size=0.8,
        validation_size=0.2,
        early_stopping=True,
        check_val_every_n_epoch=1,
        accelerator="cpu",
    )
    hist_elbo = nichevae.history["elbo_train"]
    latent = nichevae.get_latent_representation()
    assert latent.shape == (adata.n_obs, nichevae.module.n_latent)
    nichevae.save("test_nichevi", save_anndata=True, overwrite=True)
    model2 = nichevae.load("test_nichevi")
    np.testing.assert_array_equal(model2.history_["elbo_train"], hist_elbo)
    latent2 = model2.get_latent_representation()
    assert np.allclose(latent, latent2, atol=1e-5)

    nichevae.get_elbo(indices=nichevae.validation_indices)
    nichevae.get_composition_error(return_mean=False, indices=nichevae.validation_indices)
    nichevae.get_niche_error(return_mean=False, indices=nichevae.validation_indices)
    nichevae.get_normalized_expression()
    predicted_alpha = nichevae.predict_neighborhood()  # specific to nicheSCVI
    assert predicted_alpha.shape == (adata.n_obs, nichevae.n_labels)
    assert np.allclose(predicted_alpha.sum(), adata.n_obs, atol=1e-5)
    predicted_eta = nichevae.predict_niche_activation()  # specific to nicheSCVI
    assert predicted_eta.shape == (adata.n_obs, nichevae.n_labels, N_LATENT)


def test_nichevi_differential(adata):
    nicheSCVI.preprocessing_anndata(
        adata,
        k_nn=K_NN,
        **setup_kwargs,
    )

    nicheSCVI.setup_anndata(
        adata,
        layer="counts",
        batch_key="batch",
        **setup_kwargs,
    )
    nichevae = nicheSCVI(
        adata,
        prior_mixture=False,
        semisupervised=True,
        linear_classifier=True,
    )

    nichevae.train(
        max_epochs=N_EPOCHS_NICHEVI,
        train_size=0.8,
        validation_size=0.2,
        early_stopping=True,
        check_val_every_n_epoch=1,
        accelerator="cpu",
    )

    nichevae.differential_expression(
        groupby="labels",
        group1="label_1",
        group2="label_2",
        batch_correction=False,
        niche_mode=False,
        fdr_target=1,
        delta=0.5,
    )

    nichevae.differential_expression(
        groupby="labels",
        group1="label_1",
        # group2="label_2",
        batch_correction=False,
        radius=None,
        k_nn=5,
        fdr_target=1,
        delta=0.5,
    )
    nichevae.differential_expression(
        groupby="labels",
        group1="label_1",
        group2="label_2",
        batch_correction=False,
        radius=50,
        k_nn=None,
        fdr_target=[1, 1, 1, 1],
        delta=[0.5, 0.5, 0.5, 0.5],
    )
    nichevae.differential_expression(
        groupby="labels",
        group1="label_1",
        group2="label_2",
        batch_correction=False,
        radius=50,
        k_nn=None,
        fdr_target=[1, 1, 1, 1],
        delta=[0.5, 0.5, 0.5, 0.5],
    )
