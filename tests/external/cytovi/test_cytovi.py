import os

import anndata as ad
import numpy as np
import pandas as pd
import pytest

from scvi.criticism import PosteriorPredictiveCheck
from scvi.criticism._constants import METRIC_CV_GENE
from scvi.data import synthetic_iid
from scvi.external import cytovi
from scvi.external.cytovi import _model as cytovi_model
from scvi.external.cytovi._utils import log_median

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


def test_cytovi_subsample_multi_modality(adata):
    # smoke test mimicking the DiagVI spatial proteomics tutorial, which combines multiple
    # modalities into one AnnData (via `ad.concat(..., label="modality")`) and then uses
    # `cytovi.subsample(..., groupby="modality")` to balance the subsample across them before
    # scib-metrics benchmarking.
    adata_rna = adata.copy()
    adata_protein = adata.copy()

    adata_combined = ad.concat(
        [adata_rna, adata_protein],
        axis=0,
        join="inner",
        label="modality",
        keys=["rna", "protein"],
    )
    adata_combined.obs_names_make_unique()

    n_obs_group = 50
    adata_combined_sub = cytovi.subsample(
        adata_combined, n_obs=2 * n_obs_group, groupby="modality"
    )

    assert adata_combined_sub.n_obs == 2 * n_obs_group
    counts = adata_combined_sub.obs["modality"].value_counts()
    assert set(counts.index) == {"rna", "protein"}
    assert (counts == n_obs_group).all()


@pytest.mark.parametrize("protein_likelihood", ["normal", "beta"])
def test_cytovi_overlapping_protein_likelihood(overlapping_adatas, protein_likelihood):
    # regression test for https://github.com/scverse/scvi-tools/issues/4006
    # merging overlapping panels replaces missing markers with 0 in the scaled layer,
    # which is outside the support of the Beta likelihood and previously produced NaN
    # training losses.
    adata1, adata2 = overlapping_adatas
    cytovi.transform_arcsinh(adata1)
    cytovi.scale(adata1)
    cytovi.transform_arcsinh(adata2)
    cytovi.scale(adata2)
    adata_merged = cytovi.merge_batches([adata1, adata2])
    assert NAN_LAYER_KEY in adata_merged.layers

    cytovi.CYTOVI.setup_anndata(
        adata_merged,
        layer=SCALED_LAYER_KEY,
        batch_key=BATCH_KEY,
    )

    model = cytovi.CYTOVI(adata_merged, protein_likelihood=protein_likelihood)
    model.train(max_epochs=N_EPOCHS)
    assert model.is_trained
    assert np.isfinite(model.history_["elbo_train"].to_numpy(dtype=float)).all()

    imp_exp = model.get_normalized_expression()
    assert imp_exp.shape == adata_merged.shape


@pytest.mark.optional
def test_cytovi_plotting(adata):
    cytovi.plot_biaxial(
        adata, layer_key=RAW_LAYER_KEY, marker_x=adata.var_names[0], show_plot=False
    )
    cytovi.plot_histogram(adata, layer_key=RAW_LAYER_KEY, show_plot=False)


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

    assert model.posterior_predictive_sample().shape == (adata.n_obs, adata.n_vars)
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


@pytest.mark.parametrize("protein_likelihood", ["normal", "beta"])
@pytest.mark.parametrize("latent_distribution", ["normal", "ln"])
def test_cytovi_likelihood_and_latent_distribution(adata, protein_likelihood, latent_distribution):
    cytovi.transform_arcsinh(adata)
    cytovi.scale(adata)

    cytovi.CYTOVI.setup_anndata(
        adata,
        layer=SCALED_LAYER_KEY,
        batch_key=BATCH_KEY,
        sample_key=SAMPLE_KEY,
    )

    model = cytovi.CYTOVI(
        adata,
        protein_likelihood=protein_likelihood,
        latent_distribution=latent_distribution,
    )
    model.train(max_epochs=N_EPOCHS)
    assert model.is_trained
    assert np.isfinite(model.history_["elbo_train"].to_numpy(dtype=float)).all()

    latent = model.get_latent_representation()
    assert latent.shape[0] == adata.n_obs

    imp_exp = model.get_normalized_expression()
    assert imp_exp.shape == adata.shape


@pytest.mark.optional
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


@pytest.fixture
def small_model():
    # untrained model on a tiny AnnData, for unit tests that only exercise plumbing
    adata = ad.AnnData(np.random.default_rng(0).uniform(0.1, 0.9, (12, 4)).astype(np.float32))
    adata.obs["sample"] = ["a"] * 6 + ["b"] * 6
    cytovi.CYTOVI.setup_anndata(adata, sample_key="sample")
    model = cytovi.CYTOVI(adata, n_latent=2)
    model.is_trained_ = True
    return model


@pytest.mark.parametrize("lfc_clipping", [False, True])
def test_cytovi_de_forwards_test_mode(small_model, monkeypatch, lfc_clipping):
    # https://github.com/scverse/scvi-tools/pull/4037
    received = {}
    monkeypatch.setattr(cytovi_model, "_de_core", lambda *args, **kwargs: received.update(kwargs))
    small_model.differential_expression(
        groupby="sample", test_mode="three", lfc_clipping=lfc_clipping, balance_samples=False
    )
    assert received["test_mode"] == "three"


def test_cytovi_log_median():
    # https://github.com/scverse/scvi-tools/pull/4040
    # even counts average the central values in probability space: log((1 + 9) / 2)
    np.testing.assert_allclose(log_median(np.log([[1.0, 9.0]])), np.log([5.0]))
    # extreme log probabilities would under/overflow if exponentiated
    np.testing.assert_allclose(log_median([[-900.0, -900.0], [880.0, 880.0]]), [-900.0, 880.0])


def test_cytovi_aggregated_posterior_scale_is_std(small_model, monkeypatch):
    # https://github.com/scverse/scvi-tools/pull/4034
    means = np.zeros((3, 2), dtype=np.float32)
    variances = np.array([[0.25, 4.0], [9.0, 16.0], [1.0, 0.0625]], dtype=np.float32)
    # get_latent_representation(return_dist=True) returns means and variances
    monkeypatch.setattr(
        small_model, "get_latent_representation", lambda **kwargs: (means, variances)
    )
    posterior = small_model.get_aggregated_posterior(indices=np.arange(3))
    np.testing.assert_allclose(
        posterior.component_distribution.scale.cpu().numpy(), np.sqrt(variances).T
    )


def test_cytovi_uses_supplied_adata(small_model):
    # https://github.com/scverse/scvi-tools/pull/4035
    adata = small_model.adata[[9, 1, 7, 3]].copy()
    posterior = small_model.get_aggregated_posterior(adata)
    assert posterior.component_distribution.loc.shape == (2, adata.n_obs)
    log_probs = small_model.get_sample_logprobs(adata)
    assert log_probs.shape == (adata.n_obs, 2)


def test_cytovi_da_alignment(small_model, monkeypatch):
    # https://github.com/scverse/scvi-tools/pull/4036
    adata = small_model.adata
    adata.obs["condition"] = adata.obs["sample"].map({"a": "case", "b": "control"}).astype(str)

    da = small_model.differential_abundance(groupby="condition")
    assert da.index.tolist() == adata.obs_names.tolist()
    assert set(da.columns) == {"DA_case", "DA_control"}

    # grouping by the sample key itself
    da = small_model.differential_abundance(groupby="sample")
    assert set(da.columns) == {"DA_a", "DA_b"}

    # conditions are matched to score columns by sample name, not by position
    scores = pd.DataFrame({"b": -3.0, "a": -1.0}, index=adata.obs_names)
    monkeypatch.setattr(small_model, "get_sample_logprobs", lambda *args, **kwargs: scores)
    da = small_model.differential_abundance(groupby="condition")
    np.testing.assert_allclose(da["DA_case"], 2.0)

    adata.obs["condition"] = ["x", "y"] * 6
    with pytest.raises(ValueError, match="exactly one condition per sample"):
        small_model.differential_abundance(groupby="condition")
