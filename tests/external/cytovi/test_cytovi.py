import pytest

from scvi.data import synthetic_iid
from scvi.external import cytovi


from anndata import AnnData
import numpy as np



RAW_LAYER_KEY = 'raw'
SCALED_LAYER_KEY = 'scaled'
NAN_LAYER_KEY = '_nan_mask'
RNA_LAYER_KEY = 'rna'
LATENT_REP_KEY = 'X_CytoVI'
BATCH_KEY = 'batch'
LABELS_KEY = 'labels'
SAMPLE_KEY = 'sample_key'
N_EPOCHS = 2

@pytest.fixture(scope="session")
def adata():
    adata = synthetic_iid(batch_size=256,
                       n_genes=30,
                       n_proteins=0,
                       n_regions=0,
                       n_batches=2,
                       n_labels=10)

    adata.layers[RAW_LAYER_KEY] = adata.X.copy()
    adata.obs[SAMPLE_KEY] = np.random.choice(['group_a', 'group_b'], size=adata.shape[0])
    return adata

@pytest.fixture(scope="session")
def overlapping_adatas():
    adata1 = synthetic_iid(batch_size=256,
                       n_genes=30,
                       n_proteins=0,
                       n_regions=0,
                       n_batches=1,
                       n_labels=10)

    adata2 = synthetic_iid(batch_size=256,
                       n_genes=20,
                       n_proteins=0,
                       n_regions=0,
                       n_batches=1,
                       n_labels=10)

    adata1.layers[RAW_LAYER_KEY] = adata1.X.copy()
    adata2.layers[RAW_LAYER_KEY] = adata2.X.copy()

    adata1.obs_names = 'adata1_' + adata1.obs_names
    adata2.obs_names = 'adata2_' + adata2.obs_names

    return adata1, adata2

def test_cytovi_preprocess(adata, overlapping_adatas):
    cytovi.logp(adata)
    cytovi.arcsinh(adata)
    cytovi.scale(adata)
    adata_sub = cytovi.subsample(adata, n_obs=100)
    assert adata_sub.n_obs == 100

    adata1, adata2 = overlapping_adatas
    cytovi.arcsinh(adata1)
    cytovi.scale(adata1)
    cytovi.arcsinh(adata2)
    cytovi.scale(adata2)
    adata_merged = cytovi.merge_batches([adata1, adata2])
    assert NAN_LAYER_KEY in adata_merged.layers

def test_cytovi_plotting(adata):
    cytovi.biaxial(adata, layer_key=RAW_LAYER_KEY, marker_x=adata.var_names[0])
    cytovi.histogram(adata, layer_key=RAW_LAYER_KEY)


def test_cytovi(adata):
    cytovi.arcsinh(adata)
    cytovi.scale(adata)

    cytovi.CYTOVI.setup_anndata(adata,
                                layer=SCALED_LAYER_KEY,
                                batch_key=BATCH_KEY,
                                sample_key=SAMPLE_KEY,
                                )

    model = cytovi.CYTOVI(adata)

    model.train(max_epochs= N_EPOCHS)
    assert model.is_trained

    latent = model.get_latent_representation()
    assert latent.shape[0] == adata.n_obs

    imp_exp = model.get_normalized_expression()
    assert imp_exp.shape == adata.shape

    model.posterior_predictive_sample()
    da_res = model.differential_abundance()
    assert da_res.shape == (adata.n_obs, adata.obs[SAMPLE_KEY].nunique())

    model.differential_expression(groupby = SAMPLE_KEY)

    # test label informed prior
    cytovi.CYTOVI.setup_anndata(adata,
                                layer=SCALED_LAYER_KEY,
                                batch_key=BATCH_KEY,
                                sample_key=SAMPLE_KEY,
                                labels_key=LABELS_KEY
                                )

    model = cytovi.CYTOVI(adata)
    model.train(max_epochs= N_EPOCHS)




def test_cytovi_overlapping(overlapping_adatas):
    adata1, adata2 = overlapping_adatas
    cytovi.arcsinh(adata1)
    cytovi.scale(adata1)
    cytovi.arcsinh(adata2)
    cytovi.scale(adata2)
    adata_merged = cytovi.merge_batches([adata1, adata2])

    cytovi.CYTOVI.setup_anndata(adata_merged,
                                layer=SCALED_LAYER_KEY,
                                batch_key=BATCH_KEY,
                                )

    model = cytovi.CYTOVI(adata_merged)

    model.train(max_epochs= N_EPOCHS)
    assert model.is_trained

    imp_exp = model.get_normalized_expression()
    assert imp_exp.shape == adata_merged.shape

    # test label imputation
    del adata1.obs[LABELS_KEY]

    adata2.obsm[LATENT_REP_KEY] = np.random.randint(0, 1, (adata2.shape[0], model.module.n_latent))
    model_query = cytovi.CYTOVI.load_query_data(adata1, model)
    model_query.is_trained = True
    imp_cats = model_query.impute_categories_from_reference(adata2, cat_key = LABELS_KEY)
    assert imp_cats.shape[0] == adata1.n_obs

    # test RNA imputation
    adata2.layers[RNA_LAYER_KEY] = np.random.randint(0, 100, size=adata2.shape)
    adata_imp_rna = model.impute_rna_from_reference(reference_batch = '1', adata_rna = adata2, layer_key=RNA_LAYER_KEY)
    assert adata_imp_rna.shape == (adata_merged.shape[0], adata2.n_vars)

