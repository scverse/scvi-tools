import tempfile

import numpy as np
import pytest
from anndata import AnnData
from sklearn.gaussian_process import GaussianProcessClassifier

from scvi.data import synthetic_iid
from scvi.external import SCVIVA
from scvi.external.scviva.differential_expression import DifferentialExpressionResults

N_LATENT_INTRINSIC = 20
N_LATENT = 10
K_NN = 5
N_EPOCHS_SCVIVA = 2

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

    adata.obsm["qz1_m"] = np.random.normal(size=(adata.shape[0], N_LATENT_INTRINSIC))
    adata.layers["counts"] = adata.X.copy()

    return adata


def test_scviva_train(adata: AnnData):
    SCVIVA.preprocessing_anndata(
        adata,
        k_nn=K_NN,
        **setup_kwargs,
    )

    SCVIVA.setup_anndata(
        adata,
        layer="counts",
        batch_key="batch",
        **setup_kwargs,
    )
    nichevae = SCVIVA(
        adata,
        prior_mixture=False,
        semisupervised=True,
        linear_classifier=True,
    )

    nichevae.train(
        max_epochs=N_EPOCHS_SCVIVA,
        train_size=0.8,
        validation_size=0.2,
        early_stopping=True,
        check_val_every_n_epoch=1,
        accelerator="cpu",
    )

    assert nichevae.is_trained


def test_scviva_save_load(adata):
    SCVIVA.preprocessing_anndata(
        adata,
        k_nn=K_NN,
        **setup_kwargs,
    )

    SCVIVA.setup_anndata(
        adata,
        layer="counts",
        batch_key="batch",
        **setup_kwargs,
    )
    nichevae = SCVIVA(
        adata,
        prior_mixture=False,
        semisupervised=True,
        linear_classifier=True,
    )

    nichevae.train(
        max_epochs=N_EPOCHS_SCVIVA,
        train_size=0.8,
        validation_size=0.2,
        early_stopping=True,
        check_val_every_n_epoch=1,
        accelerator="cpu",
    )
    hist_elbo = nichevae.history["elbo_train"]
    latent = nichevae.get_latent_representation()
    assert latent.shape == (adata.n_obs, nichevae.module.n_latent)
    nichevae.save("test_scVIVA", save_anndata=True, overwrite=True)
    model2 = nichevae.load("test_scVIVA")
    np.testing.assert_array_equal(model2.history_["elbo_train"], hist_elbo)
    latent2 = model2.get_latent_representation()
    assert np.allclose(latent, latent2, atol=1e-5)

    nichevae.get_elbo(indices=nichevae.validation_indices)
    nichevae.get_composition_error(return_mean=False, indices=nichevae.validation_indices)
    nichevae.get_niche_error(return_mean=False, indices=nichevae.validation_indices)
    nichevae.get_normalized_expression()
    predicted_alpha = nichevae.predict_neighborhood()
    assert predicted_alpha.shape == (adata.n_obs, nichevae.n_labels)
    assert np.allclose(predicted_alpha.sum(), adata.n_obs, atol=1e-5)
    predicted_eta = nichevae.predict_niche_activation()
    assert predicted_eta.shape == (adata.n_obs, nichevae.n_labels, N_LATENT_INTRINSIC)


def test_scviva_differential(adata):
    SCVIVA.preprocessing_anndata(
        adata,
        k_nn=K_NN,
        **setup_kwargs,
    )

    SCVIVA.setup_anndata(
        adata,
        layer="counts",
        batch_key="batch",
        **setup_kwargs,
    )
    nichevae = SCVIVA(
        adata,
        prior_mixture=False,
        semisupervised=True,
        linear_classifier=True,
    )

    nichevae.train(
        max_epochs=N_EPOCHS_SCVIVA,
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
    DE_results = nichevae.differential_expression(
        groupby="labels",
        group1="label_1",
        group2="label_2",
        batch_correction=False,
        radius=None,
        k_nn=5,
        fdr_target=[1, 1, 1, 1],
        delta=[0.5, 0.5, 0.5, 0.5],
    )

    assert isinstance(DE_results, DifferentialExpressionResults)
    assert isinstance(DE_results.gpc, GaussianProcessClassifier)
    assert hasattr(DE_results.gpc, "log_marginal_likelihood_value_")

    # Suppress plt.show() to avoid UI popups
    import matplotlib.pyplot as plt

    plt_show_backup = plt.show
    plt.show = lambda: None

    try:
        DE_results.plot()
    finally:
        plt.show = plt_show_backup

    import os

    with tempfile.NamedTemporaryFile(suffix=".png", delete=False) as tmp:
        path = tmp.name

    try:
        DE_results.plot(path_to_save=path)
        assert os.path.exists(path)
        assert os.path.getsize(path) > 0
    finally:
        os.remove(path)

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


def test_scviva_scarches(adata: AnnData):
    # divide between ref and query data
    adata.obs["hemisphere"] = [
        "right" if x > 0 else "left" for x in adata.obsm["coordinates"][:, 0]
    ]
    ref_adata = adata[adata.obs["hemisphere"] == "left"].copy()
    query_adata = adata[adata.obs["hemisphere"] == "right"].copy()

    SCVIVA.preprocessing_anndata(
        ref_adata,
        k_nn=K_NN,
        **setup_kwargs,
    )

    SCVIVA.setup_anndata(
        ref_adata,
        layer="counts",
        batch_key="batch",
        **setup_kwargs,
    )

    # Reference model
    nichevae = SCVIVA(
        ref_adata,
        prior_mixture=False,
        semisupervised=True,
        linear_classifier=True,
    )

    nichevae.train(
        max_epochs=N_EPOCHS_SCVIVA,
        train_size=0.8,
        validation_size=0.2,
        early_stopping=True,
        check_val_every_n_epoch=1,
        accelerator="cpu",
    )

    assert nichevae.is_trained

    # Query adata
    SCVIVA.preprocessing_anndata(
        query_adata,
        k_nn=K_NN,
        **setup_kwargs,
    )

    # Make it different from reference adata - How?
    query_adata = query_adata[
        :, np.random.permutation(query_adata.var_names)[: query_adata.shape[1] - 5]
    ].copy()
    query_adata.obsm["neighborhood_composition"] = query_adata.obsm[
        "neighborhood_composition"
    ].reindex(columns=ref_adata.obsm["neighborhood_composition"].columns)

    # Make query adata and model
    nichevae.prepare_query_anndata(query_adata, reference_model=nichevae)
    query_nichevae = nichevae.load_query_data(query_adata, reference_model=nichevae)

    query_nichevae.train(
        max_epochs=N_EPOCHS_SCVIVA,
        train_size=0.8,
        validation_size=0.2,
        early_stopping=True,
        check_val_every_n_epoch=1,
        accelerator="cpu",
    )

    predicted_alpha = query_nichevae.predict_neighborhood(query_adata)
    assert predicted_alpha.shape == (query_adata.n_obs, query_nichevae.n_labels)
    query_adata.obsm["X_nichevi"] = query_nichevae.get_latent_representation(query_adata)

    # Check that embedding default works
    assert (
        query_nichevae.get_latent_representation(
            adata=ref_adata,
        ).shape[0]
        == ref_adata.shape[0]
    )
    assert (
        query_nichevae.get_latent_representation(
            adata=query_adata,
        ).shape[0]
        == query_adata.shape[0]
    )
