import os
import random

import anndata
import numpy as np
import pandas as pd
import pytest
from scipy import sparse
from scipy.sparse.csr import csr_matrix

import scvi
from scvi import REGISTRY_KEYS
from scvi.data import _constants, synthetic_iid
from scvi.data.fields import ObsmField, ProteinObsmField
from scvi.dataloaders import AnnTorchDataset

from .utils import generic_setup_adata_manager


def test_transfer_fields_basic(adata1, adata2):
    # test transfer_fields function
    adata2.X = adata1.X
    adata1_manager = generic_setup_adata_manager(adata1)
    adata1_manager.transfer_fields(adata2)
    np.testing.assert_array_equal(
        adata1.obs["_scvi_labels"], adata2.obs["_scvi_labels"]
    )


def test_transfer_fields_layer_use(adata1, adata2):
    # test if layer was used initially, again used in transfer setup
    zeros = np.zeros_like(adata1.X)
    ones = np.ones_like(adata1.X)
    adata1.X = zeros
    adata2.X = ones
    adata1_manager = generic_setup_adata_manager(adata1, layer="raw")
    adata1_manager.transfer_fields(adata2)
    np.testing.assert_array_equal(
        adata1.obs["_scvi_labels"], adata2.obs["_scvi_labels"]
    )


def test_transfer_fields_unknown_batch(adata1, adata2):
    # test that an unknown batch throws an error
    adata2.obs["batch"] = [2] * adata2.n_obs
    adata1_manager = generic_setup_adata_manager(adata1, batch_key="batch")
    with pytest.raises(ValueError):
        adata1_manager.transfer_fields(adata2)


def test_transfer_fields_unknown_label(adata1, adata2):
    # test that an unknown label throws an error
    adata2.obs["labels"] = ["label_123"] * adata2.n_obs
    adata1_manager = generic_setup_adata_manager(adata1, labels_key="labels")
    with pytest.raises(ValueError):
        adata1_manager.transfer_fields(adata2)


def test_transfer_fields_correct_mapping(adata1, adata2):
    # test that correct mapping was applied
    adata2.obs["labels"] = ["label_1"] * adata2.n_obs
    adata1_manager = generic_setup_adata_manager(adata1, labels_key="labels")
    adata1_manager.transfer_fields(adata2)
    labels_mapping = adata1_manager.get_state_registry("labels").categorical_mapping
    correct_label = np.where(labels_mapping == "label_1")[0][0]
    assert adata2.obs["_scvi_labels"][0] == correct_label


def test_transfer_fields_correct_batch(adata1, adata2):
    # test that transfer_fields correctly looks for adata.obs['batch']
    del adata2.obs["batch"]
    adata1_manager = generic_setup_adata_manager(adata1, batch_key="batch")
    with pytest.raises(KeyError):
        adata1_manager.transfer_fields(adata2)


def test_transfer_fields_same_batch_and_label(adata1, adata2):
    # test that transfer_fields assigns same batch and label to cells
    # if the original anndata was also same batch and label
    adata1_manager = generic_setup_adata_manager(adata1)
    del adata2.obs["batch"]
    adata1_manager.transfer_fields(adata2)
    assert adata2.obs["_scvi_batch"][0] == 0
    assert adata2.obs["_scvi_labels"][0] == 0


def test_transfer_fields_subset(adata1, adata2):
    # test that if a category mapping is a subset, transfer anndata is called
    scvi.model.SCVI.setup_anndata(adata1, batch_key="batch")
    adata2.obs["batch"] = "batch_1"
    scvi.model.SCVI.setup_anndata(adata2, batch_key="batch")
    m = scvi.model.SCVI(adata1)
    m.train(1)
    m.get_latent_representation(adata2)
    assert adata2.obs["_scvi_batch"].all() == 1


def test_transfer_fields_wrong_kwarg(adata):
    # test that error is thrown if an arbitrary kwarg is passed into setup_anndata
    with pytest.raises(TypeError):
        scvi.model.SCVI.setup_anndata(adata, batch="batch")


def test_clobber_same_model(adata):
    scvi.model.SCVI.setup_anndata(adata)
    m1 = scvi.model.SCVI(adata)
    m1.train(1)

    scvi.model.SCVI.setup_anndata(adata, batch_key="batch")
    m2 = scvi.model.SCVI(adata)
    m2.train(1)

    adata_manager1 = m1.get_anndata_manager(adata)
    assert adata_manager1.summary_stats.n_batch == 1
    # The underlying data is still 2 since we have not run _validate_anndata yet
    # to re-transfer the setup of m1.
    assert (
        len(np.unique(adata_manager1.get_from_registry(REGISTRY_KEYS.BATCH_KEY))) == 2
    )
    m1._validate_anndata(adata)
    assert (
        len(np.unique(adata_manager1.get_from_registry(REGISTRY_KEYS.BATCH_KEY))) == 1
    )

    adata_manager2 = m2.get_anndata_manager(adata)
    assert adata_manager2.summary_stats.n_batch == 2
    # The underlying data is still 1 since we have not run _validate_anndata yet
    # to re-transfer the setup of m2.
    assert (
        len(np.unique(adata_manager2.get_from_registry(REGISTRY_KEYS.BATCH_KEY))) == 1
    )
    m2._validate_anndata(adata)
    assert (
        len(np.unique(adata_manager2.get_from_registry(REGISTRY_KEYS.BATCH_KEY))) == 2
    )


def test_clobber_different_models(adata):
    scvi.model.SCVI.setup_anndata(adata, batch_key="batch")
    m1 = scvi.model.SCVI(adata)
    m1.train(1)

    scvi.model.TOTALVI.setup_anndata(
        adata,
        protein_expression_obsm_key="protein_expression",
        protein_names_uns_key="protein_names",
    )
    m2 = scvi.model.TOTALVI(adata)
    m2.train(1)

    adata_manager1 = m1.get_anndata_manager(adata)
    assert adata_manager1.summary_stats.n_batch == 2
    # The underlying data is still 2 since we have not run _validate_anndata yet
    # to re-transfer the setup of m1.
    assert (
        len(np.unique(adata_manager1.get_from_registry(REGISTRY_KEYS.BATCH_KEY))) == 1
    )
    m1._validate_anndata(adata)
    assert (
        len(np.unique(adata_manager1.get_from_registry(REGISTRY_KEYS.BATCH_KEY))) == 2
    )

    adata_manager2 = m2.get_anndata_manager(adata)
    assert adata_manager2.summary_stats.n_batch == 1
    # The underlying data is still 1 since we have not run _validate_anndata yet
    # to re-transfer the setup of m2.
    assert (
        len(np.unique(adata_manager2.get_from_registry(REGISTRY_KEYS.BATCH_KEY))) == 2
    )
    m2._validate_anndata(adata)
    assert (
        len(np.unique(adata_manager2.get_from_registry(REGISTRY_KEYS.BATCH_KEY))) == 1
    )


def test_register_new_fields(adata):
    bdata = adata.copy()
    incremental_adata_manager = generic_setup_adata_manager(bdata)
    batch_field = None
    for field in incremental_adata_manager.fields:
        if field.registry_key == REGISTRY_KEYS.BATCH_KEY:
            batch_field = field
    new_fields = [
        ProteinObsmField(
            REGISTRY_KEYS.PROTEIN_EXP_KEY,
            "protein_expression",
            batch_field=batch_field,
            use_batch_mask=True,
            is_count_data=True,
        )
    ]
    incremental_adata_manager.register_new_fields(new_fields)
    adata_manager = generic_setup_adata_manager(
        adata, protein_expression_obsm_key="protein_expression"
    )
    np.testing.assert_array_equal(
        incremental_adata_manager.get_from_registry(REGISTRY_KEYS.PROTEIN_EXP_KEY),
        adata_manager.get_from_registry(REGISTRY_KEYS.PROTEIN_EXP_KEY),
    )
    assert len(incremental_adata_manager.fields) == len(adata_manager.fields)


def test_register_new_fields_with_transferred_manager(adata):
    bdata = adata.copy()
    cdata = adata.copy()
    adata_manager = generic_setup_adata_manager(adata, batch_key="batch")
    bdata.obs["batch"] = adata.obs["batch"].to_numpy()[0]
    bdata_manager = adata_manager.transfer_fields(bdata)
    new_fields = [ObsmField(REGISTRY_KEYS.PROTEIN_EXP_KEY, "protein_expression")]
    bdata_manager.register_new_fields(new_fields)
    cdata_manager = bdata_manager.transfer_fields(cdata)

    # Should have protein field
    cdata_manager.get_from_registry(REGISTRY_KEYS.PROTEIN_EXP_KEY)
    np.testing.assert_array_equal(
        cdata.obs["_scvi_batch"].values, adata.obs["_scvi_batch"].values
    )


def test_update_setup_args(adata):
    adata_manager = generic_setup_adata_manager(adata)
    adata_manager.update_setup_method_args({"test_arg": "test_val"})
    assert (
        "test_arg"
        in adata_manager._get_setup_method_args()[_constants._SETUP_ARGS_KEY].keys()
    )


def test_data_format(adata):
    # if data was dense np array, check after setup_anndata, data is C_CONTIGUOUS
    old_x = adata.X
    old_pro = adata.obsm["protein_expression"]
    old_obs = adata.obs
    adata.X = np.asfortranarray(old_x)
    adata.obsm["protein_expression"] = np.asfortranarray(old_pro)
    assert adata.X.flags["C_CONTIGUOUS"] is False
    assert adata.obsm["protein_expression"].flags["C_CONTIGUOUS"] is False

    adata_manager = generic_setup_adata_manager(
        adata, protein_expression_obsm_key="protein_expression"
    )
    assert adata.X.flags["C_CONTIGUOUS"] is True
    assert adata.obsm["protein_expression"].flags["C_CONTIGUOUS"] is True

    assert np.array_equal(old_x, adata.X)
    assert np.array_equal(old_pro, adata.obsm["protein_expression"])
    assert np.array_equal(old_obs, adata.obs)

    assert np.array_equal(adata.X, adata_manager.get_from_registry(REGISTRY_KEYS.X_KEY))
    assert np.array_equal(
        adata.obsm["protein_expression"],
        adata_manager.get_from_registry(REGISTRY_KEYS.PROTEIN_EXP_KEY),
    )


def test_data_format_c_contiguous(adata):
    # if obsm is dataframe, make it C_CONTIGUOUS if it isnt
    pe = np.asfortranarray(adata.obsm["protein_expression"])
    adata.obsm["protein_expression"] = pd.DataFrame(pe, index=adata.obs_names)
    assert adata.obsm["protein_expression"].to_numpy().flags["C_CONTIGUOUS"] is False
    adata_manager = generic_setup_adata_manager(
        adata, protein_expression_obsm_key="protein_expression"
    )
    new_pe = adata_manager.get_from_registry(REGISTRY_KEYS.PROTEIN_EXP_KEY)
    assert new_pe.to_numpy().flags["C_CONTIGUOUS"] is True
    assert np.array_equal(pe, new_pe)
    assert np.array_equal(adata.X, adata_manager.get_from_registry(REGISTRY_KEYS.X_KEY))
    assert np.array_equal(
        adata.obsm["protein_expression"],
        adata_manager.get_from_registry(REGISTRY_KEYS.PROTEIN_EXP_KEY),
    )


def test_setup_anndata(adata):
    # test regular setup
    adata_manager = generic_setup_adata_manager(
        adata,
        batch_key="batch",
        labels_key="labels",
        protein_expression_obsm_key="protein_expression",
        protein_names_uns_key="protein_names",
    )
    np.testing.assert_array_equal(
        adata_manager.get_from_registry(REGISTRY_KEYS.BATCH_KEY),
        np.array(adata.obs["_scvi_batch"]).reshape((-1, 1)),
    )
    np.testing.assert_array_equal(
        adata_manager.get_from_registry(REGISTRY_KEYS.LABELS_KEY),
        np.array(adata.obs["labels"].cat.codes).reshape((-1, 1)),
    )
    np.testing.assert_array_equal(
        adata_manager.get_from_registry(REGISTRY_KEYS.X_KEY), adata.X
    )
    np.testing.assert_array_equal(
        adata_manager.get_from_registry(REGISTRY_KEYS.PROTEIN_EXP_KEY),
        adata.obsm["protein_expression"],
    )
    np.testing.assert_array_equal(
        adata_manager.get_state_registry(REGISTRY_KEYS.PROTEIN_EXP_KEY).column_names,
        adata.uns["protein_names"],
    )


def test_setup_anndata_view_error(adata):
    # test that error is thrown if its a view:
    with pytest.raises(ValueError):
        generic_setup_adata_manager(adata[1])


def test_setup_anndata_view_error_df_protein_none(adata):
    # If obsm is a df and protein_names_uns_key is None, protein names should be grabbed from column of df
    new_protein_names = np.array(random.sample(range(100), 100)).astype("str")
    df = pd.DataFrame(
        adata.obsm["protein_expression"],
        index=adata.obs_names,
        columns=new_protein_names,
    )
    adata.obsm["protein_expression"] = df
    adata_manager = generic_setup_adata_manager(
        adata, protein_expression_obsm_key="protein_expression"
    )
    np.testing.assert_array_equal(
        adata_manager.get_state_registry(REGISTRY_KEYS.PROTEIN_EXP_KEY).column_names,
        new_protein_names,
    )


def test_setup_anndata_layer(adata):
    # test that layer is working properly
    true_x = adata.X
    adata.layers["X"] = true_x
    adata.X = np.ones_like(adata.X)
    adata_manager = generic_setup_adata_manager(adata, layer="X")
    np.testing.assert_array_equal(
        adata_manager.get_from_registry(REGISTRY_KEYS.X_KEY), true_x
    )


def test_setup_anndata_create_label_batch(adata):
    # test that it creates labels and batch if no layers_key is passed
    adata_manager = generic_setup_adata_manager(
        adata,
        protein_expression_obsm_key="protein_expression",
        protein_names_uns_key="protein_names",
    )
    np.testing.assert_array_equal(
        adata_manager.get_from_registry(REGISTRY_KEYS.BATCH_KEY),
        np.zeros((adata.shape[0], 1)),
    )
    np.testing.assert_array_equal(
        adata_manager.get_from_registry(REGISTRY_KEYS.LABELS_KEY),
        np.zeros((adata.shape[0], 1)),
    )


def test_setup_anndata_nan(adata):
    # test error is thrown when categorical obs field contains nans
    adata.obs["batch"][:10] = np.nan
    with pytest.raises(ValueError):
        generic_setup_adata_manager(adata, batch_key="batch")


def test_setup_anndata_cat(adata):
    # test error is thrown when categorical joint obsm field contains nans
    adata.obs["cat1"][:10] = np.nan
    with pytest.raises(ValueError):
        generic_setup_adata_manager(adata, categorical_covariate_keys=["cat1"])


def test_save_setup_anndata(adata, save_path):
    generic_setup_adata_manager(
        adata,
        batch_key="batch",
        labels_key="labels",
        protein_expression_obsm_key="protein_expression",
        protein_names_uns_key="protein_names",
    )
    adata.write(os.path.join(save_path, "test.h5ad"))


def test_extra_covariates(adata):
    scvi.model.SCVI.setup_anndata(
        adata,
        batch_key="batch",
        labels_key="labels",
        continuous_covariate_keys=["cont1", "cont2"],
        categorical_covariate_keys=["cat1", "cat2"],
    )
    m = scvi.model.SCVI(adata)
    m.train(1)
    df1 = m.get_from_registry(adata, REGISTRY_KEYS.CONT_COVS_KEY)
    df2 = adata.obs[["cont1", "cont2"]]
    pd.testing.assert_frame_equal(df1, df2)


def test_extra_covariates_transfer(adata):
    adata_manager = generic_setup_adata_manager(
        adata,
        batch_key="batch",
        labels_key="labels",
        protein_expression_obsm_key="protein_expression",
        protein_names_uns_key="protein_names",
        continuous_covariate_keys=["cont1", "cont2"],
        categorical_covariate_keys=["cat1", "cat2"],
    )
    bdata = synthetic_iid()
    bdata.obs["cont1"] = np.random.normal(size=(bdata.shape[0],))
    bdata.obs["cont2"] = np.random.normal(size=(bdata.shape[0],))
    bdata.obs["cat1"] = 0
    bdata.obs["cat2"] = 1

    adata_manager.transfer_fields(bdata)

    # give it a new category
    bdata.obs["cat1"] = 6
    bdata_manager = adata_manager.transfer_fields(bdata, extend_categories=True)
    assert (
        bdata_manager.get_state_registry(REGISTRY_KEYS.CAT_COVS_KEY).mappings["cat1"][
            -1
        ]
        == 6
    )


def test_anntorchdataset_getitem(adata):
    adata_manager = generic_setup_adata_manager(
        adata,
        batch_key="batch",
        labels_key="labels",
        protein_expression_obsm_key="protein_expression",
        protein_names_uns_key="protein_names",
    )
    # check that we can successfully pass in a list of tensors to get
    tensors_to_get = [REGISTRY_KEYS.BATCH_KEY, REGISTRY_KEYS.LABELS_KEY]
    bd = AnnTorchDataset(adata_manager, getitem_tensors=tensors_to_get)
    np.testing.assert_array_equal(tensors_to_get, list(bd[1].keys()))

    # check that we can successfully pass in a dict of tensors and their associated types
    bd = AnnTorchDataset(
        adata_manager,
        getitem_tensors={
            REGISTRY_KEYS.X_KEY: np.int,
            REGISTRY_KEYS.LABELS_KEY: np.int64,
        },
    )
    assert bd[1][REGISTRY_KEYS.X_KEY].dtype == np.int64
    assert bd[1][REGISTRY_KEYS.LABELS_KEY].dtype == np.int64

    # check that by default we get all the registered tensors
    bd = AnnTorchDataset(adata_manager)
    all_registered_tensors = list(adata_manager.data_registry.keys())
    np.testing.assert_array_equal(all_registered_tensors, list(bd[1].keys()))
    assert bd[1][REGISTRY_KEYS.X_KEY].shape[0] == bd.adata_manager.summary_stats.n_vars


def test_anntorchdataset_numpy(adata):
    # check that AnnTorchDataset returns numpy array
    adata_manager = generic_setup_adata_manager(adata)
    bd = AnnTorchDataset(adata_manager)
    for value in bd[1].values():
        assert type(value) == np.ndarray


def test_anntorchdataset_numpy_sparse(adata):
    # check AnnTorchDataset returns numpy array counts were sparse
    adata.X = sparse.csr_matrix(adata.X)
    adata_manager = generic_setup_adata_manager(adata)
    bd = AnnTorchDataset(adata_manager)
    for value in bd[1].values():
        assert type(value) == np.ndarray


def test_anntorchdataset_getitem_numpy_sparse(adata):
    # check AnnTorchDataset returns numpy array if pro exp was sparse
    adata.obsm["protein_expression"] = sparse.csr_matrix(
        adata.obsm["protein_expression"]
    )
    adata_manager = generic_setup_adata_manager(
        adata, batch_key="batch", protein_expression_obsm_key="protein_expression"
    )
    bd = AnnTorchDataset(adata_manager)
    for value in bd[1].values():
        assert type(value) == np.ndarray


def test_anntorchdataset_getitem_pro_exp(adata):
    # check pro exp is being returned as numpy array even if its DF
    adata.obsm["protein_expression"] = pd.DataFrame(
        adata.obsm["protein_expression"], index=adata.obs_names
    )
    adata_manager = generic_setup_adata_manager(
        adata, batch_key="batch", protein_expression_obsm_key="protein_expression"
    )
    bd = AnnTorchDataset(adata_manager)
    for value in bd[1].values():
        assert type(value) == np.ndarray


def test_view_registry(adata):
    adata_manager = generic_setup_adata_manager(
        adata,
        batch_key="batch",
        labels_key="labels",
        protein_expression_obsm_key="protein_expression",
        protein_names_uns_key="protein_names",
        continuous_covariate_keys=["cont1", "cont2"],
        categorical_covariate_keys=["cat1", "cat2"],
    )
    adata_manager.view_registry()
    adata_manager.view_registry(hide_state_registries=True)


def test_saving(adata, save_path):
    save_path = os.path.join(save_path, "tmp_adata.h5ad")
    adata.obs["cat1"] = np.random.randint(0, 3, adata.n_obs).astype(str)
    adata.obs["cat1"][1] = "asdf"
    adata.obs["cat1"][2] = "f34"
    adata.obs["cat2"] = np.random.randint(0, 7, adata.n_obs)

    generic_setup_adata_manager(
        adata,
        protein_expression_obsm_key="protein_expression",
        batch_key="batch",
        labels_key="labels",
        categorical_covariate_keys=["cat1", "cat2"],
        continuous_covariate_keys=["cont1", "cont2"],
    )
    adata.write(save_path)
    anndata.read(save_path)


def test_backed_anndata(adata, save_path):
    path = os.path.join(save_path, "test_data.h5ad")
    adata.write_h5ad(path)
    adata = anndata.read_h5ad(path, backed="r+")
    adata_manager = generic_setup_adata_manager(adata, batch_key="batch")

    # test get item
    bd = AnnTorchDataset(adata_manager)
    subset = bd[np.arange(adata.n_obs)]
    assert isinstance(subset["X"], np.ndarray)


def test_backed_anndata_sparse(adata, save_path):
    # sparse
    adata.X = csr_matrix(adata.X)
    path = os.path.join(save_path, "test_data2.h5ad")
    adata.write_h5ad(path)
    adata = anndata.read_h5ad(path, backed="r+")
    adata_manager = generic_setup_adata_manager(adata, batch_key="batch")

    # test get item
    bd = AnnTorchDataset(adata_manager)
    subset = bd[np.arange(adata.n_obs)]
    assert isinstance(subset["X"], np.ndarray)
