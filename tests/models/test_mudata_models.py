import os

import numpy as np
import pytest
from anndata import AnnData
from mudata import MuData

import scvi
from scvi.data import synthetic_iid
from scvi.model import TOTALVI
from scvi.utils import attrdict


def test_totalvi(save_path):
    adata = synthetic_iid()
    protein_adata = synthetic_iid(n_genes=50)
    mdata = MuData({"rna": adata, "protein": protein_adata})
    TOTALVI.setup_mudata(
        mdata,
        batch_key="batch",
        modalities=dict(rna_layer="rna", batch_key="rna", protein_layer="protein"),
    )

    n_obs = mdata.n_obs
    n_genes = adata.n_vars
    n_proteins = protein_adata.X.shape[1]
    n_latent = 10

    model = TOTALVI(mdata, n_latent=n_latent)
    model.train(1, train_size=0.5)
    assert model.is_trained is True
    z = model.get_latent_representation()
    assert z.shape == (n_obs, n_latent)
    model.get_elbo()
    model.get_marginal_ll(n_mc_samples=3)
    model.get_reconstruction_error()
    model.get_normalized_expression()
    model.get_normalized_expression(transform_batch=["batch_0", "batch_1"])
    model.get_latent_library_size()
    model.get_protein_foreground_probability()
    model.get_protein_foreground_probability(transform_batch=["batch_0", "batch_1"])
    post_pred = model.posterior_predictive_sample(n_samples=2)
    assert post_pred.shape == (n_obs, n_genes + n_proteins, 2)
    post_pred = model.posterior_predictive_sample(n_samples=1)
    assert post_pred.shape == (n_obs, n_genes + n_proteins)
    feature_correlation_matrix1 = model.get_feature_correlation_matrix(
        correlation_type="spearman"
    )
    feature_correlation_matrix1 = model.get_feature_correlation_matrix(
        correlation_type="spearman", transform_batch=["batch_0", "batch_1"]
    )
    feature_correlation_matrix2 = model.get_feature_correlation_matrix(
        correlation_type="pearson"
    )
    assert feature_correlation_matrix1.shape == (
        n_genes + n_proteins,
        n_genes + n_proteins,
    )
    assert feature_correlation_matrix2.shape == (
        n_genes + n_proteins,
        n_genes + n_proteins,
    )

    model.get_elbo(indices=model.validation_indices)
    model.get_marginal_ll(indices=model.validation_indices, n_mc_samples=3)
    model.get_reconstruction_error(indices=model.validation_indices)

    adata2 = synthetic_iid()
    protein_adata2 = synthetic_iid(n_genes=50)
    mdata2 = MuData({"rna": adata, "protein": protein_adata})
    TOTALVI.setup_mudata(
        mdata2,
        batch_key="batch",
        modalities=dict(rna_layer="rna", batch_key="rna", protein_layer="protein"),
    )
    norm_exp = model.get_normalized_expression(mdata2, indices=[1, 2, 3])
    assert norm_exp[0].shape == (3, adata2.n_vars)
    assert norm_exp[1].shape == (3, protein_adata2.n_vars)
    norm_exp = model.get_normalized_expression(
        mdata2,
        gene_list=adata2.var_names[:5].to_list(),
        protein_list=protein_adata2.var_names[:3].to_list(),
        transform_batch=["batch_0", "batch_1"],
    )

    latent_lib_size = model.get_latent_library_size(mdata2, indices=[1, 2, 3])
    assert latent_lib_size.shape == (3, 1)

    pro_foreground_prob = model.get_protein_foreground_probability(
        mdata2, indices=[1, 2, 3], protein_list=["1", "2"]
    )
    assert pro_foreground_prob.shape == (3, 2)
    model.posterior_predictive_sample(mdata2)
    model.get_feature_correlation_matrix(mdata2)

    # test transfer_anndata_setup + view
    adata2 = synthetic_iid()
    protein_adata2 = synthetic_iid(n_genes=50)
    mdata2 = MuData({"rna": adata2, "protein": protein_adata2})
    model.get_elbo(mdata2[:10])


def test_totalvi_auto_transfer(save_path):
    # test automatic transfer_fields
    adata = synthetic_iid()
    protein_adata = synthetic_iid(n_genes=50)
    mdata = MuData({"rna": adata, "protein": protein_adata})
    TOTALVI.setup_mudata(
        mdata,
        batch_key="batch",
        modalities=dict(rna_layer="rna", batch_key="rna", protein_layer="protein"),
    )
    model = TOTALVI(mdata)
    adata2 = synthetic_iid()
    protein_adata2 = synthetic_iid(n_genes=50)
    mdata2 = MuData({"rna": adata2, "protein": protein_adata2})
    model.get_elbo(mdata2)


def test_totalvi_incorrect_mapping(save_path):
    # test that we catch incorrect mappings
    adata = synthetic_iid()
    protein_adata = synthetic_iid(n_genes=50)
    mdata = MuData({"rna": adata, "protein": protein_adata})
    TOTALVI.setup_mudata(
        mdata,
        batch_key="batch",
        modalities=dict(rna_layer="rna", batch_key="rna", protein_layer="protein"),
    )
    model = TOTALVI(mdata)
    adata2 = synthetic_iid()
    protein_adata2 = synthetic_iid(n_genes=50)
    mdata2 = MuData({"rna": adata2, "protein": protein_adata2})
    adata2.obs.batch.cat.rename_categories(["batch_0", "batch_10"], inplace=True)
    with pytest.raises(ValueError):
        model.get_elbo(mdata2)


def test_totalvi_reordered_mapping(save_path):
    # test that same mapping different order is okay
    adata = synthetic_iid()
    protein_adata = synthetic_iid(n_genes=50)
    mdata = MuData({"rna": adata, "protein": protein_adata})
    TOTALVI.setup_mudata(
        mdata,
        batch_key="batch",
        modalities=dict(rna_layer="rna", batch_key="rna", protein_layer="protein"),
    )
    model = TOTALVI(mdata)
    adata2 = synthetic_iid()
    protein_adata2 = synthetic_iid(n_genes=50)
    mdata2 = MuData({"rna": adata2, "protein": protein_adata2})
    adata2.obs.batch.cat.rename_categories(["batch_1", "batch_0"], inplace=True)
    model.get_elbo(mdata2)


def test_totalvi_missing_proteins(save_path):
    # test with missing proteins
    adata = scvi.data.pbmcs_10x_cite_seq(
        save_path=save_path,
        protein_join="outer",
    )
    protein_adata = AnnData(adata.obsm["protein_expression"])
    mdata = MuData({"rna": adata, "protein": protein_adata})
    TOTALVI.setup_mudata(
        mdata,
        batch_key="batch",
        modalities=dict(rna_layer="rna", batch_key="rna", protein_layer="protein"),
    )
    model = TOTALVI(mdata)
    assert model.module.protein_batch_mask is not None
    model.train(1, train_size=0.5)

    model = TOTALVI(mdata, override_missing_proteins=True)
    assert model.module.protein_batch_mask is None
    model.train(1, train_size=0.5)


def test_totalvi_model_library_size(save_path):
    adata = synthetic_iid()
    protein_adata = synthetic_iid(n_genes=50)
    mdata = MuData({"rna": adata, "protein": protein_adata})
    TOTALVI.setup_mudata(
        mdata,
        batch_key="batch",
        modalities=dict(rna_layer="rna", batch_key="rna", protein_layer="protein"),
    )

    n_latent = 10
    model = TOTALVI(mdata, n_latent=n_latent, use_observed_lib_size=False)
    assert hasattr(model.module, "library_log_means") and hasattr(
        model.module, "library_log_vars"
    )
    model.train(1, train_size=0.5)
    assert model.is_trained is True
    model.get_elbo()
    model.get_marginal_ll(n_mc_samples=3)
    model.get_latent_library_size()


def test_totalvi_size_factor():
    adata = synthetic_iid()
    adata.obs["size_factor"] = np.random.randint(1, 5, size=(adata.shape[0],))
    protein_adata = synthetic_iid(n_genes=50)
    mdata = MuData({"rna": adata, "protein": protein_adata})
    TOTALVI.setup_mudata(
        mdata,
        batch_key="batch",
        size_factor_key="size_factor",
        modalities=dict(
            rna_layer="rna",
            batch_key="rna",
            protein_layer="protein",
            size_factor_key="rna",
        ),
    )

    n_latent = 10

    # Test size_factor_key overrides use_observed_lib_size.
    model = TOTALVI(mdata, n_latent=n_latent, use_observed_lib_size=False)
    assert not hasattr(model.module, "library_log_means") and not hasattr(
        model.module, "library_log_vars"
    )
    assert model.module.use_size_factor_key
    model.train(1, train_size=0.5)

    model = TOTALVI(mdata, n_latent=n_latent, use_observed_lib_size=True)
    assert not hasattr(model.module, "library_log_means") and not hasattr(
        model.module, "library_log_vars"
    )
    assert model.module.use_size_factor_key
    model.train(1, train_size=0.5)


def test_totalvi_saving_and_loading(save_path):
    adata = synthetic_iid()
    protein_adata = synthetic_iid(n_genes=50)
    mdata = MuData({"rna": adata, "protein": protein_adata})
    TOTALVI.setup_mudata(
        mdata,
        batch_key="batch",
        modalities=dict(rna_layer="rna", batch_key="rna", protein_layer="protein"),
    )
    model = TOTALVI(mdata)
    model.train(1, train_size=0.2)
    z1 = model.get_latent_representation(mdata)
    test_idx1 = model.validation_indices

    model.save(save_path, overwrite=True, save_anndata=True)
    model.view_setup_args(save_path)

    model = TOTALVI.load(save_path)
    model.get_latent_representation()

    # Load with mismatched genes.
    tmp_adata = synthetic_iid(
        n_genes=200,
    )
    tmp_protein_adata = synthetic_iid(n_genes=50)
    tmp_mdata = MuData({"rna": tmp_adata, "protein": tmp_protein_adata})
    with pytest.raises(ValueError):
        TOTALVI.load(save_path, adata=tmp_mdata)

    # Load with different batches.
    tmp_adata = synthetic_iid()
    tmp_adata.obs["batch"] = tmp_adata.obs["batch"].cat.rename_categories(
        ["batch_2", "batch_3"]
    )
    tmp_protein_adata = synthetic_iid(n_genes=50)
    tmp_mdata = MuData({"rna": tmp_adata, "protein": tmp_protein_adata})
    with pytest.raises(ValueError):
        TOTALVI.load(save_path, adata=tmp_mdata)

    model = TOTALVI.load(save_path, adata=mdata)
    assert scvi.REGISTRY_KEYS.BATCH_KEY in model.adata_manager.data_registry
    assert model.adata_manager.data_registry.batch == attrdict(
        dict(mod_key="rna", attr_name="obs", attr_key="_scvi_batch")
    )

    z2 = model.get_latent_representation()
    test_idx2 = model.validation_indices
    np.testing.assert_array_equal(z1, z2)
    np.testing.assert_array_equal(test_idx1, test_idx2)
    assert model.is_trained is True

    save_path = os.path.join(save_path, "tmp")

    adata2 = synthetic_iid()
    protein_adata2 = synthetic_iid(n_genes=50)
    mdata2 = MuData({"rna": adata2, "protein": protein_adata2})
    TOTALVI.setup_mudata(
        mdata2,
        batch_key="batch",
        modalities=dict(rna_layer="rna", batch_key="rna", protein_layer="protein"),
    )
