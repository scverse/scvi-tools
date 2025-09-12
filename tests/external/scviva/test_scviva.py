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
LABELS_KEY = "labels"

setup_kwargs = {
    "sample_key": "batch",
    "labels_key": LABELS_KEY,
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
        n_batches=3,
        n_labels=5,
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


@pytest.fixture
def split_ref_query_adata(adata: AnnData):
    """Split adata into ref and query using hemisphere."""
    adata.obs["hemisphere"] = [
        "right" if x > 0 else "left" for x in adata.obsm["coordinates"][:, 0]
    ]
    ref_adata = adata[adata.obs["hemisphere"] == "left"].copy()
    query_adata = adata[adata.obs["hemisphere"] == "right"].copy()
    return ref_adata, query_adata


def test_scviva_scarches_less_features(split_ref_query_adata):
    # divide between ref and query data
    ref_adata, query_adata = split_ref_query_adata

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

    # Make it different from reference adata - shuffling and removing a few genes here
    query_adata = query_adata[
        :, np.random.permutation(query_adata.var_names)[: query_adata.shape[1] - 5]
    ].copy()

    query_adata.obs["labels"] = (
        query_adata.obs["labels"].astype(str).replace("label_2", "label_1").astype("category")
    )
    # Query adata
    nichevae.preprocessing_query_anndata(
        query_adata,
        reference_model=nichevae,
        k_nn=K_NN,
        **setup_kwargs,
    )

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


def _permute_label_categories(adata: AnnData, label_key: str = "labels") -> None:
    """Permute the order of label categories to ensure .unique() returns different order."""
    import random

    current_categories = list(adata.obs[label_key].cat.categories)
    new_order = random.sample(current_categories, len(current_categories))  # shuffled list
    adata.obs[label_key] = adata.obs[label_key].cat.reorder_categories(new_order, ordered=False)


def test_scviva_scarches_same_features(split_ref_query_adata):
    # divide between ref and query data
    ref_adata, query_adata = split_ref_query_adata

    ref_labels = set(ref_adata.obs[LABELS_KEY].unique())
    query_labels = set(query_adata.obs[LABELS_KEY].unique())
    assert ref_labels == query_labels, (
        f"Label sets do not match:\nRef: {ref_labels}\nQuery: {query_labels}"
    )

    ref_adata.obs["labels"] = ref_adata.obs[LABELS_KEY].astype("category")
    _permute_label_categories(ref_adata, label_key=LABELS_KEY)
    query_adata.obs["labels"] = query_adata.obs[LABELS_KEY].astype("category")
    _permute_label_categories(query_adata, label_key=LABELS_KEY)

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

    # Make it different from reference adata - only shuffling genes here
    query_adata = query_adata[
        :, np.random.permutation(query_adata.var_names)[: query_adata.shape[1]]
    ].copy()

    # Query adata
    nichevae.preprocessing_query_anndata(
        query_adata,
        reference_model=nichevae,
        k_nn=K_NN,
        **setup_kwargs,
    )

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
