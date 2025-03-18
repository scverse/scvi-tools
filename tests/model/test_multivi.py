import os

import anndata as ad
import numpy as np
import pytest
import scanpy as sc
from mudata import MuData

import scvi
from scvi import REGISTRY_KEYS
from scvi.data import synthetic_iid
from scvi.model import MULTIVI
from scvi.utils import attrdict


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
    vae.get_normalized_accessibility()
    vae.get_normalized_accessibility(normalize_cells=True)
    vae.get_normalized_accessibility(normalize_regions=True)
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
    MULTIVI.setup_anndata(data, batch_key="batch")
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
    assert vae.n_proteins == data.obsm["protein_expression"].shape[1]
    vae.train(3)


def test_multivi_single_batch():
    data = synthetic_iid(n_batches=1)
    MULTIVI.setup_anndata(data, batch_key="batch")
    vae = MULTIVI(data, n_genes=50, n_regions=50)
    with pytest.warns(UserWarning):
        vae.train(3)


@pytest.mark.internet
def test_multivi_mudata_rna_prot_external():
    # Example on how to download protein adata to mudata (from multivi tutorial) - mudata RNA/PROT
    adata = scvi.data.pbmcs_10x_cite_seq()
    adata.layers["counts"] = adata.X.copy()
    adata.obs_names_make_unique()
    protein_adata = ad.AnnData(adata.obsm["protein_expression"])
    protein_adata.obs_names = adata.obs_names
    del adata.obsm["protein_expression"]
    mdata = MuData({"rna": adata, "protein": protein_adata})
    sc.pp.highly_variable_genes(
        mdata.mod["rna"],
        n_top_genes=4000,
        flavor="seurat_v3",
        batch_key="batch",
        layer="counts",
    )
    mdata.mod["rna_subset"] = mdata.mod["rna"][:, mdata.mod["rna"].var["highly_variable"]].copy()
    mdata.update()
    # mdata
    # mdata.mod
    MULTIVI.setup_mudata(
        mdata,
        rna_layer="counts",  # mean we use: mdata.mod["rna_subset"].layers["counts"]
        protein_layer=None,  # mean we use: mdata.mod["protein"].X
        batch_key="batch",  # the batch is here: mdata.mod["rna_subset"].obs["batch"]
        modalities={
            "rna_layer": "rna_subset",
            "protein_layer": "protein",
            "batch_key": "rna_subset",
        },
    )
    model = MULTIVI(mdata)
    model.train(1, train_size=0.9)


def test_multivi_mudata_rna_atac_external():
    # optional data - mudata RNA/ATAC
    mdata = synthetic_iid(return_mudata=True)
    sc.pp.highly_variable_genes(
        mdata.mod["rna"],
        n_top_genes=4000,
        flavor="seurat_v3",
    )
    mdata.mod["rna_subset"] = mdata.mod["rna"][:, mdata.mod["rna"].var["highly_variable"]].copy()
    sc.pp.highly_variable_genes(
        mdata.mod["accessibility"],
        n_top_genes=4000,
        flavor="seurat_v3",
    )
    mdata.mod["atac_subset"] = mdata.mod["accessibility"][
        :, mdata.mod["accessibility"].var["highly_variable"]
    ].copy()
    mdata.update()
    MULTIVI.setup_mudata(
        mdata,
        modalities={"rna_layer": "rna_subset", "atac_layer": "atac_subset"},
    )
    model = MULTIVI(mdata)
    model.train(1, train_size=0.9)


def test_multivi_mudata_trimodal_external():
    # optional data - mudata RNA/ATAC
    mdata = synthetic_iid(return_mudata=True)
    MULTIVI.setup_mudata(
        mdata,
        modalities={
            "rna_layer": "rna",
            "atac_layer": "accessibility",
            "protein_layer": "protein_expression",
        },
    )
    model = MULTIVI(mdata)
    model.train(1, train_size=0.9)
    model.train(1, train_size=0.9)
    assert model.is_trained is True
    model.get_latent_representation()
    model.get_elbo()
    model.get_reconstruction_error()
    model.get_normalized_expression()
    model.get_normalized_accessibility()
    model.get_normalized_accessibility(normalize_cells=True)
    model.get_normalized_accessibility(normalize_regions=True)
    model.get_library_size_factors()
    model.get_region_factors()

    model.get_elbo(indices=model.validation_indices)
    model.get_reconstruction_error(indices=model.validation_indices)
    model.get_normalized_accessibility()
    model.get_normalized_accessibility(normalize_cells=True)
    model.get_normalized_accessibility(normalize_regions=True)
    model.get_library_size_factors()
    model.get_region_factors()


@pytest.mark.parametrize("n_genes", [25, 50, 100])
@pytest.mark.parametrize("n_regions", [25, 50, 100])
def test_multivi_mudata(n_genes: int, n_regions: int):
    # use of syntetic data of rna/proteins/atac for speed

    mdata = synthetic_iid(return_mudata=True)
    MULTIVI.setup_mudata(
        mdata,
        batch_key="batch",
        modalities={
            "rna_layer": "rna",
            "protein_layer": "protein_expression",
            "atac_layer": "accessibility",
        },
    )
    n_obs = mdata.n_obs
    n_latent = 10

    model = MULTIVI(mdata, n_latent=n_latent, n_genes=n_genes, n_regions=n_regions)
    model.train(1, train_size=0.9)
    assert model.is_trained is True
    z = model.get_latent_representation()
    assert z.shape == (n_obs, n_latent)
    model.get_elbo()
    model.get_reconstruction_error()
    model.get_normalized_expression()
    model.get_normalized_expression(transform_batch=["batch_0", "batch_1"])
    model.get_normalized_accessibility()
    model.get_normalized_accessibility(normalize_cells=True)
    model.get_normalized_accessibility(normalize_regions=True)
    model.get_library_size_factors()
    model.get_region_factors()

    model.get_elbo(indices=model.validation_indices)
    model.get_reconstruction_error(indices=model.validation_indices)
    model.get_normalized_accessibility()
    model.get_normalized_accessibility(normalize_cells=True)
    model.get_normalized_accessibility(normalize_regions=True)
    model.get_library_size_factors()
    model.get_region_factors()

    mdata2 = synthetic_iid(return_mudata=True)
    MULTIVI.setup_mudata(
        mdata2,
        batch_key="batch",
        modalities={"rna_layer": "rna", "protein_layer": "protein_expression"},
    )
    norm_exp = model.get_normalized_expression(mdata2, indices=[1, 2, 3])
    assert norm_exp.shape == (3, n_genes)

    # test transfer_anndata_setup + view
    mdata3 = synthetic_iid(return_mudata=True)
    mdata3.obs["_indices"] = np.arange(mdata3.n_obs)
    model.get_elbo(mdata3[:10])
    model.get_normalized_accessibility()
    model.get_normalized_accessibility(normalize_cells=True)
    model.get_normalized_accessibility(normalize_regions=True)
    model.get_library_size_factors()
    model.get_region_factors()


def test_multivi_auto_transfer_mudata():
    # test automatic transfer_fields
    adata = synthetic_iid()
    protein_adata = synthetic_iid(n_genes=50)
    mdata = MuData({"rna": adata, "protein": protein_adata})
    MULTIVI.setup_mudata(
        mdata,
        batch_key="batch",
        modalities={"rna_layer": "rna", "batch_key": "rna", "protein_layer": "protein"},
    )
    model = MULTIVI(mdata)
    adata2 = synthetic_iid()
    protein_adata2 = synthetic_iid(n_genes=50)
    mdata2 = MuData({"rna": adata2, "protein": protein_adata2})
    mdata2.obs["_indices"] = np.arange(mdata2.n_obs)
    model.get_elbo(mdata2)
    model.get_normalized_accessibility()
    model.get_normalized_accessibility(normalize_cells=True)
    model.get_normalized_accessibility(normalize_regions=True)
    model.get_library_size_factors()
    model.get_region_factors()


def test_multivi_incorrect_mapping_mudata():
    # test that we catch incorrect mappings
    adata = synthetic_iid()
    protein_adata = synthetic_iid(n_genes=50)
    mdata = MuData({"rna": adata, "protein": protein_adata})
    MULTIVI.setup_mudata(
        mdata,
        batch_key="batch",
        modalities={"rna_layer": "rna", "batch_key": "rna", "protein_layer": "protein"},
    )
    model = MULTIVI(mdata)
    adata2 = synthetic_iid()
    protein_adata2 = synthetic_iid(n_genes=50)
    mdata2 = MuData({"rna": adata2, "protein": protein_adata2})
    adata2.obs.batch = adata2.obs.batch.cat.rename_categories(["batch_0", "batch_10"])
    with pytest.raises(ValueError):
        model.get_elbo(mdata2)


def test_multivi_reordered_mapping_mudata():
    # test that same mapping different order is okay
    adata = synthetic_iid()
    protein_adata = synthetic_iid(n_genes=50)
    mdata = MuData({"rna": adata, "protein": protein_adata})
    MULTIVI.setup_mudata(
        mdata,
        batch_key="batch",
        modalities={"rna_layer": "rna", "batch_key": "rna", "protein_layer": "protein"},
    )
    model = MULTIVI(mdata)
    adata2 = synthetic_iid()
    protein_adata2 = synthetic_iid(n_genes=50)
    mdata2 = MuData({"rna": adata2, "protein": protein_adata2})
    adata2.obs.batch = adata2.obs.batch.cat.rename_categories(["batch_1", "batch_0"])
    mdata2.obs["_indices"] = np.arange(mdata2.n_obs)
    model.get_elbo(mdata2)
    model.get_normalized_accessibility()
    model.get_normalized_accessibility(normalize_cells=True)
    model.get_normalized_accessibility(normalize_regions=True)
    model.get_library_size_factors()
    model.get_region_factors()


def test_multivi_model_library_size_mudata():
    adata = synthetic_iid()
    protein_adata = synthetic_iid(n_genes=50)
    mdata = MuData({"rna": adata, "protein": protein_adata})
    MULTIVI.setup_mudata(
        mdata,
        batch_key="batch",
        modalities={"rna_layer": "rna", "batch_key": "rna", "protein_layer": "protein"},
    )

    n_latent = 10
    model = MULTIVI(mdata, n_latent=n_latent)
    model.train(1, train_size=0.5)
    assert model.is_trained is True
    model.get_elbo()
    model.get_normalized_accessibility()
    model.get_normalized_accessibility(normalize_cells=True)
    model.get_normalized_accessibility(normalize_regions=True)
    model.get_library_size_factors()
    model.get_region_factors()


def test_multivi_size_factor_mudata():
    mdata = synthetic_iid(return_mudata=True)
    mdata.obs["size_factor_rna"] = mdata["rna"].X.sum(1)
    mdata.obs["size_factor_atac"] = (mdata["accessibility"].X.sum(1) + 1) / (
        np.max(mdata["accessibility"].X.sum(1)) + 1.01
    )
    MULTIVI.setup_mudata(
        mdata,
        modalities={"rna_layer": "rna", "atac_layer": "accessibility"},
        size_factor_key=["size_factor_rna", "size_factor_atac"],
    )

    n_latent = 10

    # Test size_factor_key overrides use_observed_lib_size.
    model = MULTIVI(mdata, n_latent=n_latent)
    assert model.module.use_size_factor_key
    model.train(1, train_size=0.5)

    model = MULTIVI(mdata, n_latent=n_latent)
    assert model.module.use_size_factor_key
    model.train(1, train_size=0.5)


def test_multivi_saving_and_loading_mudata(save_path: str):
    adata = synthetic_iid()
    protein_adata = synthetic_iid(n_genes=50)
    mdata = MuData({"rna": adata, "protein": protein_adata})
    MULTIVI.setup_mudata(
        mdata,
        batch_key="batch",
        modalities={"rna_layer": "rna", "batch_key": "rna", "protein_layer": "protein"},
    )
    model = MULTIVI(mdata)
    model.train(1, train_size=0.2)
    z1 = model.get_latent_representation(mdata)
    test_idx1 = model.validation_indices

    model.save(save_path, overwrite=True, save_anndata=True)
    model.view_setup_args(save_path)

    model = MULTIVI.load(save_path)
    model.get_latent_representation()

    # Load with mismatched genes.
    tmp_adata = synthetic_iid(
        n_genes=200,
    )
    tmp_protein_adata = synthetic_iid(n_genes=50)
    tmp_mdata = MuData({"rna": tmp_adata, "protein": tmp_protein_adata})
    with pytest.raises(ValueError):
        MULTIVI.load(save_path, adata=tmp_mdata)

    # Load with different batches.
    tmp_adata = synthetic_iid()
    tmp_adata.obs["batch"] = tmp_adata.obs["batch"].cat.rename_categories(["batch_2", "batch_3"])
    tmp_protein_adata = synthetic_iid(n_genes=50)
    tmp_mdata = MuData({"rna": tmp_adata, "protein": tmp_protein_adata})
    with pytest.raises(ValueError):
        MULTIVI.load(save_path, adata=tmp_mdata)

    model = MULTIVI.load(save_path, adata=mdata)
    assert REGISTRY_KEYS.BATCH_KEY in model.adata_manager.data_registry
    assert model.adata_manager.data_registry.batch == attrdict(
        {"mod_key": "rna", "attr_name": "obs", "attr_key": "_scvi_batch"}
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
    MULTIVI.setup_mudata(
        mdata2,
        batch_key="batch",
        modalities={"rna_layer": "rna", "batch_key": "rna", "protein_layer": "protein"},
    )


def test_scarches_mudata_prep_layer(save_path: str):
    n_latent = 5
    mdata1 = synthetic_iid(return_mudata=True)

    mdata1["rna"].layers["counts"] = mdata1["rna"].X.copy()
    MULTIVI.setup_mudata(
        mdata1,
        batch_key="batch",
        modalities={"rna_layer": "rna", "protein_layer": "protein_expression"},
    )
    model = MULTIVI(mdata1, n_latent=n_latent)
    model.train(1, check_val_every_n_epoch=1)
    dir_path = os.path.join(save_path, "saved_model/")
    model.save(dir_path, overwrite=True)

    # mdata2 has more genes and missing 10 genes from mdata1.
    # protein/acessibility features are same as in mdata1
    mdata2 = synthetic_iid(n_genes=110, return_mudata=True)
    mdata2["rna"].layers["counts"] = mdata2["rna"].X.copy()
    new_var_names_init = [f"Random {i}" for i in range(10)]
    new_var_names = new_var_names_init + mdata2["rna"].var_names[10:].to_list()
    mdata2["rna"].var_names = new_var_names

    original_protein_values = mdata2["protein_expression"].X.copy()
    original_accessibility_values = mdata2["accessibility"].X.copy()

    MULTIVI.prepare_query_mudata(mdata2, dir_path)
    # should be padded 0s
    assert np.sum(mdata2["rna"][:, mdata2["rna"].var_names[:10]].layers["counts"]) == 0
    np.testing.assert_equal(
        mdata2["rna"].var_names[:10].to_numpy(), mdata1["rna"].var_names[:10].to_numpy()
    )

    # values of other modalities should be unchanged
    np.testing.assert_equal(original_protein_values, mdata2["protein_expression"].X)
    np.testing.assert_equal(original_accessibility_values, mdata2["accessibility"].X)

    # and names should also be the same
    np.testing.assert_equal(
        mdata2["protein_expression"].var_names.to_numpy(),
        mdata1["protein_expression"].var_names.to_numpy(),
    )
    np.testing.assert_equal(
        mdata2["accessibility"].var_names.to_numpy(), mdata1["accessibility"].var_names.to_numpy()
    )
    MULTIVI.load_query_data(mdata2, dir_path)


def test_multivi_save_load_mudata_format(save_path: str):
    mdata = synthetic_iid(return_mudata=True, protein_expression_key="protein")
    invalid_mdata = mdata.copy()
    invalid_mdata.mod["protein"] = invalid_mdata.mod["protein"][:, :10].copy()
    MULTIVI.setup_mudata(
        mdata,
        modalities={"rna_layer": "rna", "protein_layer": "protein"},
    )
    model = MULTIVI(mdata)
    model.train(max_epochs=1)

    legacy_model_path = os.path.join(save_path, "legacy_model")
    model.save(
        legacy_model_path,
        overwrite=True,
        save_anndata=False,
        legacy_mudata_format=True,
    )

    with pytest.raises(ValueError):
        _ = MULTIVI.load(legacy_model_path, adata=invalid_mdata)
    model = MULTIVI.load(legacy_model_path, adata=mdata)

    model_path = os.path.join(save_path, "model")
    model.save(
        model_path,
        overwrite=True,
        save_anndata=False,
        legacy_mudata_format=False,
    )
    with pytest.raises(ValueError):
        _ = MULTIVI.load(legacy_model_path, adata=invalid_mdata)
    model = MULTIVI.load(model_path, adata=mdata)
