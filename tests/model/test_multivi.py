import muon
import numpy as np
import pytest
from mudata import MuData

from scvi.data import synthetic_iid
from scvi.model import MULTIVI


def test_multivi():
    data = synthetic_iid()
    MULTIVI.setup_anndata(
        data,
        batch_key="batch",
    )
    vae = MULTIVI(
        data,
        n_genes=50,
        n_regions=50,
    )
    vae.train(1, save_best=False)
    vae.train(1, adversarial_mixing=False)
    vae.train(3)
    vae.get_elbo(indices=vae.validation_indices)
    vae.get_accessibility_estimates()
    vae.get_accessibility_estimates(normalize_cells=True)
    vae.get_accessibility_estimates(normalize_regions=True)
    vae.get_normalized_expression()
    vae.get_library_size_factors()
    vae.get_region_factors()
    vae.get_reconstruction_error(indices=vae.validation_indices)
    vae.get_latent_representation()
    vae.differential_accessibility(groupby="labels", group1="label_1")
    vae.differential_expression(groupby="labels", group1="label_1")

    # Test with size factor
    data = synthetic_iid()
    data.obs["size_factor"] = np.random.randint(1, 5, size=(data.shape[0],))
    MULTIVI.setup_anndata(data, batch_key="batch", size_factor_key="size_factor")
    vae = MULTIVI(
        data,
        n_genes=50,
        n_regions=50,
    )
    vae.train(3)

    # Test with modality weights and penalties
    data = synthetic_iid()
    MULTIVI.setup_anndata(data, batch_key="batch")
    vae = MULTIVI(data, n_genes=50, n_regions=50, modality_weights="cell")
    vae.train(3)
    vae = MULTIVI(data, n_genes=50, n_regions=50, modality_weights="universal")
    vae.train(3)
    vae = MULTIVI(data, n_genes=50, n_regions=50, modality_penalty="MMD")
    vae.train(3)

    # Test with non-zero protein data
    data = synthetic_iid()
    MULTIVI.setup_anndata(
        data,
        batch_key="batch",
        protein_expression_obsm_key="protein_expression",
        protein_names_uns_key="protein_names",
    )
    vae = MULTIVI(
        data,
        n_genes=50,
        n_regions=50,
        modality_weights="cell",
    )
    # assert vae.n_regions == data.obsm["protein_expression"].shape[1]
    vae.train(3)


def test_multivi_single_batch():
    data = synthetic_iid(n_batches=1)
    MULTIVI.setup_anndata(data, batch_key="batch")
    vae = MULTIVI(data, n_genes=50, n_regions=50)
    with pytest.warns(UserWarning):
        vae.train(3)


def test_multivi_mudata():
    # optional data - big one
    url = ("https://cf.10xgenomics.com/samples/cell-arc/2.0.0/10k_PBMC_Multiome_nextgem_Chromium_X"
           "/10k_PBMC_Multiome_nextgem_Chromium_X_filtered_feature_bc_matrix.h5")
    mdata = muon.read_10x_h5("data/multiome10k.h5mu", backup_url=url)
    mdata
    MULTIVI.setup_mudata(mdata, modalities={"rna_layer": "rna", "protein_layer": "atac"})
    vae = MULTIVI(mdata, n_genes=50, n_regions=50)

    adata = synthetic_iid()
    protein_adata = synthetic_iid(n_genes=50)
    mdata = MuData({"rna": adata, "protein": protein_adata})
    MULTIVI.setup_mudata(
        mdata,
        batch_key="batch",
        modalities={"rna_layer": "rna", "batch_key": "rna", "protein_layer": "protein"},
    )

    n_obs = mdata.n_obs
    n_genes = np.min([adata.n_vars, protein_adata.n_vars])
    n_regions = protein_adata.X.shape[1]
    n_latent = 10

    model = MULTIVI(mdata, n_latent=n_latent, n_genes=n_genes, n_regions=n_regions)
    model.train(1, train_size=0.5)
    assert model.is_trained is True
    z = model.get_latent_representation()
    assert z.shape == (n_obs, n_latent)
    model.get_elbo()
    # model.get_marginal_ll(n_mc_samples=3)
    model.get_reconstruction_error()
    model.get_normalized_expression()
    model.get_normalized_expression(transform_batch=["batch_0", "batch_1"])
    # model.get_latent_library_size()
    # model.get_protein_foreground_probability()
    # model.get_protein_foreground_probability(transform_batch=["batch_0", "batch_1"])
    # post_pred = model.posterior_predictive_sample(n_samples=2)
    # assert post_pred.shape == (n_obs, n_genes + n_regions, 2)
    # post_pred = model.posterior_predictive_sample(n_samples=1)
    # assert post_pred.shape == (n_obs, n_genes + n_regions)
    # feature_correlation_matrix1 = model.get_feature_correlation_matrix
    #                               (correlation_type="spearman")
    # feature_correlation_matrix1 = model.get_feature_correlation_matrix(
    #     correlation_type="spearman", transform_batch=["batch_0", "batch_1"]
    # )
    # feature_correlation_matrix2 = (model.get_feature_correlation_matrix
    #                                (correlation_type="pearson"))
    # assert feature_correlation_matrix1.shape == (
    #     n_genes + n_regions,
    #     n_genes + n_regions,
    # )
    # assert feature_correlation_matrix2.shape == (
    #     n_genes + n_regions,
    #     n_genes + n_regions,
    # )

    model.get_elbo(indices=model.validation_indices)
    # model.get_marginal_ll(indices=model.validation_indices, n_mc_samples=3)
    model.get_reconstruction_error(indices=model.validation_indices)

    adata2 = synthetic_iid()
    protein_adata2 = synthetic_iid(n_genes=50)
    mdata2 = MuData({"rna": adata, "protein": protein_adata})
    MULTIVI.setup_mudata(
        mdata2,
        batch_key="batch",
        modalities={"rna_layer": "rna", "batch_key": "rna", "protein_layer": "protein"},
    )
    norm_exp = model.get_normalized_expression(mdata2, indices=[1, 2, 3])
    # assert norm_exp[0].shape == (3, adata2.n_vars)
    # assert norm_exp[1].shape == (3, protein_adata2.n_vars)
    # norm_exp = model.get_normalized_expression(
    #     mdata2,
    #     gene_list=adata2.var_names[:5].to_list(),
    #     protein_list=protein_adata2.var_names[:3].to_list(),
    #     transform_batch=["batch_0", "batch_1"],
    # )

    # latent_lib_size = model.get_latent_library_size(mdata2, indices=[1, 2, 3])
    # assert latent_lib_size.shape == (3, 1)

    # pro_foreground_prob = model.get_protein_foreground_probability(
    #     mdata2, indices=[1, 2, 3], protein_list=["gene_1", "gene_2"]
    # )
    # assert pro_foreground_prob.shape == (3, 2)
    # model.posterior_predictive_sample(mdata2)
    # model.get_feature_correlation_matrix(mdata2)

    # test transfer_anndata_setup + view
    adata2 = synthetic_iid()
    protein_adata2 = synthetic_iid(n_genes=50)
    mdata2 = MuData({"rna": adata2, "protein": protein_adata2})
    # model.get_elbo(mdata2[:10])
