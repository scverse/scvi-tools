import os

import numpy as np
import pytest

from scvi.criticism import PosteriorPredictiveCheck
from scvi.criticism._constants import METRIC_CV_GENE
from scvi.data import synthetic_iid
from scvi.external import cytovi

RAW_LAYER_KEY = "raw"
SCALED_LAYER_KEY = "scaled"
NAN_LAYER_KEY = "_nan_mask"
RNA_LAYER_KEY = "rna"
LATENT_REP_KEY = "X_CytoVI"
BATCH_KEY = "batch"
LABELS_KEY = "labels"
SAMPLE_KEY = "sample_key"
N_EPOCHS = 2


@pytest.fixture(scope="session")
def adata():
    adata = synthetic_iid(
        batch_size=256,
        n_genes=30,
        n_proteins=0,
        n_regions=0,
        n_batches=2,
        n_labels=10,
        rna_dist="normal",
    )

    adata.layers[RAW_LAYER_KEY] = adata.X.copy()
    adata.obs[SAMPLE_KEY] = np.random.choice(["group_a", "group_b"], size=adata.shape[0])
    return adata


@pytest.fixture(scope="session")
def overlapping_adatas():
    adata1 = synthetic_iid(
        batch_size=256,
        n_genes=30,
        n_proteins=0,
        n_regions=0,
        n_batches=1,
        n_labels=10,
        rna_dist="normal",
    )

    adata2 = synthetic_iid(
        batch_size=256,
        n_genes=20,
        n_proteins=0,
        n_regions=0,
        n_batches=1,
        n_labels=10,
        rna_dist="normal",
    )

    adata1.layers[RAW_LAYER_KEY] = adata1.X.copy()
    adata2.layers[RAW_LAYER_KEY] = adata2.X.copy()

    adata1.obs_names = "adata1_" + adata1.obs_names
    adata2.obs_names = "adata2_" + adata2.obs_names

    return adata1, adata2


def test_cytovi_preprocess(adata, overlapping_adatas):
    cytovi.transform_arcsinh(adata)
    cytovi.scale(adata)
    adata_sub = cytovi.subsample(adata, n_obs=100)
    assert adata_sub.n_obs == 100

    adata1, adata2 = overlapping_adatas
    cytovi.transform_arcsinh(adata1)
    cytovi.scale(adata1)
    cytovi.transform_arcsinh(adata2)
    cytovi.scale(adata2)
    adata_merged = cytovi.merge_batches([adata1, adata2])
    assert NAN_LAYER_KEY in adata_merged.layers


def test_cytovi_plotting(adata):
    cytovi.plot_biaxial(adata, layer_key=RAW_LAYER_KEY, marker_x=adata.var_names[0])
    cytovi.plot_histogram(adata, layer_key=RAW_LAYER_KEY)


def test_cytovi(adata):
    cytovi.transform_arcsinh(adata)
    cytovi.scale(adata)

    cytovi.CYTOVI.setup_anndata(
        adata,
        layer=SCALED_LAYER_KEY,
        batch_key=BATCH_KEY,
        sample_key=SAMPLE_KEY,
    )

    model = cytovi.CYTOVI(adata)

    model.train(max_epochs=N_EPOCHS)
    assert model.is_trained

    ppc = PosteriorPredictiveCheck(
        adata, models_dict={"mymodel": model, "mymodel2": model}, count_layer_key="scaled"
    )
    ppc.coefficient_of_variation()
    assert ppc.metrics[METRIC_CV_GENE].shape[0] == adata.n_vars

    latent = model.get_latent_representation()
    assert latent.shape[0] == adata.n_obs

    imp_exp = model.get_normalized_expression()
    assert imp_exp.shape == adata.shape

    model.posterior_predictive_sample()
    da_res = model.differential_abundance()
    assert da_res.shape == (adata.n_obs, adata.obs[SAMPLE_KEY].nunique())

    model.differential_expression(groupby=SAMPLE_KEY)

    # test label informed prior
    cytovi.CYTOVI.setup_anndata(
        adata,
        layer=SCALED_LAYER_KEY,
        batch_key=BATCH_KEY,
        sample_key=SAMPLE_KEY,
        labels_key=LABELS_KEY,
    )

    model = cytovi.CYTOVI(adata)
    model.train(max_epochs=N_EPOCHS)


def test_cytovi_overlapping(overlapping_adatas):
    adata1, adata2 = overlapping_adatas
    cytovi.transform_arcsinh(adata1)
    cytovi.scale(adata1)
    cytovi.transform_arcsinh(adata2)
    cytovi.scale(adata2)
    adata_merged = cytovi.merge_batches([adata1, adata2])

    cytovi.CYTOVI.setup_anndata(
        adata_merged,
        layer=SCALED_LAYER_KEY,
        batch_key=BATCH_KEY,
    )

    model = cytovi.CYTOVI(adata_merged)

    model.train(max_epochs=N_EPOCHS)
    assert model.is_trained

    imp_exp = model.get_normalized_expression()
    assert imp_exp.shape == adata_merged.shape

    # test label imputation
    del adata1.obs[LABELS_KEY]

    adata2.obsm[LATENT_REP_KEY] = np.random.randint(0, 1, (adata2.shape[0], model.module.n_latent))
    model_query = cytovi.CYTOVI.load_query_data(adata1, model)
    model_query.is_trained = True
    imp_cats = model_query.impute_categories_from_reference(adata2, cat_key=LABELS_KEY)
    assert imp_cats.shape[0] == adata1.n_obs

    # test RNA imputation
    adata2.layers[RNA_LAYER_KEY] = np.random.randint(0, 100, size=adata2.shape)
    adata_imp_rna = model.impute_rna_from_reference(
        reference_batch="1", adata_rna=adata2, layer_key=RNA_LAYER_KEY
    )
    assert adata_imp_rna.shape == (adata_merged.shape[0], adata2.n_vars)


def test_cytovi_save_load(adata, save_path):
    cytovi.transform_arcsinh(adata)
    cytovi.scale(adata)

    cytovi.CYTOVI.setup_anndata(
        adata,
        layer=SCALED_LAYER_KEY,
        batch_key=BATCH_KEY,
        sample_key=SAMPLE_KEY,
    )

    model = cytovi.CYTOVI(adata)

    model.train(max_epochs=N_EPOCHS)
    hist_elbo = model.history_["elbo_train"]
    latent = model.get_latent_representation()
    assert latent.shape == (adata.n_obs, model.module.n_latent)

    model_path = os.path.join(save_path, "test_cytovi")

    model.save(model_path, save_anndata=True, overwrite=True)
    model2 = model.load(model_path)
    np.testing.assert_array_equal(model2.history_["elbo_train"], hist_elbo)
    latent2 = model2.get_latent_representation()
    assert np.allclose(latent, latent2)


def test_cytovi_write_read_fcs(adata, save_path):
    cytovi.transform_arcsinh(adata)
    cytovi.scale(adata)

    cytovi.write_fcs(adata, output_path=save_path, prefix="test_cytovi", layer=SCALED_LAYER_KEY)
    adata_read = cytovi.read_fcs(save_path + "test_cytovi.fcs")
    assert adata_read.shape == adata.shape
    assert np.allclose(adata_read.X, adata.layers[SCALED_LAYER_KEY])
