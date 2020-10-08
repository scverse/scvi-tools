import numpy as np
import os
import random
import pandas as pd
import scipy.sparse as sparse

import pytest
from scvi.data._scvidataset import ScviDataset
from scvi.data import synthetic_iid
from scvi.data import (
    setup_anndata,
    transfer_anndata_setup,
    register_tensor_from_anndata,
)
from scvi import _CONSTANTS
from scvi.data._anndata import get_from_registry
from scipy.sparse import csc_matrix


def test_transfer_anndata_setup():
    # test if raw was initially used, that it is again used in transfer
    adata1 = synthetic_iid(run_setup_anndata=False)
    adata2 = synthetic_iid(run_setup_anndata=False)
    raw_counts = adata1.copy()
    adata1.raw = raw_counts
    adata2.raw = raw_counts
    zeros = np.zeros_like(adata1.X)
    ones = np.ones_like(adata1.X)
    adata1.X = zeros
    adata2.X = ones
    setup_anndata(adata1, use_raw=True)
    transfer_anndata_setup(adata1, adata2)
    # makes sure that use_raw was used
    np.testing.assert_array_equal(
        adata1.obs["_scvi_local_l_mean"], adata2.obs["_scvi_local_l_mean"]
    )

    # test if layer was used initially, again used in transfer setup
    adata1 = synthetic_iid(run_setup_anndata=False)
    adata2 = synthetic_iid(run_setup_anndata=False)
    raw_counts = adata1.X.copy()
    adata1.layers["raw"] = raw_counts
    adata2.layers["raw"] = raw_counts
    zeros = np.zeros_like(adata1.X)
    ones = np.ones_like(adata1.X)
    adata1.X = zeros
    adata2.X = ones
    setup_anndata(adata1, layer="raw")
    transfer_anndata_setup(adata1, adata2)
    np.testing.assert_array_equal(
        adata1.obs["_scvi_local_l_mean"], adata2.obs["_scvi_local_l_mean"]
    )

    # test that an unknown batch throws an error
    adata1 = synthetic_iid()
    adata2 = synthetic_iid(run_setup_anndata=False)
    adata2.obs["batch"] = [2] * adata2.n_obs
    with pytest.raises(ValueError):
        transfer_anndata_setup(adata1, adata2)

    # TODO: test that a batch with wrong dtype throws an error
    # adata1 = synthetic_iid()
    # adata2 = synthetic_iid(run_setup_anndata=False)
    # adata2.obs["batch"] = ["0"] * adata2.n_obs
    # with pytest.raises(ValueError):
    #     transfer_anndata_setup(adata1, adata2)

    # test that an unknown label throws an error
    adata1 = synthetic_iid()
    adata2 = synthetic_iid(run_setup_anndata=False)
    adata2.obs["labels"] = ["label_123"] * adata2.n_obs
    with pytest.raises(ValueError):
        transfer_anndata_setup(adata1, adata2)

    # test that correct mapping was applied
    adata1 = synthetic_iid()
    adata2 = synthetic_iid(run_setup_anndata=False)
    adata2.obs["labels"] = ["label_1"] * adata2.n_obs
    transfer_anndata_setup(adata1, adata2)
    labels_mapping = adata1.uns["_scvi"]["categorical_mappings"]["_scvi_labels"][
        "mapping"
    ]
    correct_label = np.where(labels_mapping == "label_1")[0][0]
    adata2.obs["_scvi_labels"][0] == correct_label

    # test that transfer_anndata_setup correctly looks for adata.obs['batch']
    adata1 = synthetic_iid()
    adata2 = synthetic_iid(run_setup_anndata=False)
    del adata2.obs["batch"]
    with pytest.raises(KeyError):
        transfer_anndata_setup(adata1, adata2)

    # test that transfer_anndata_setup assigns same batch and label to cells
    # if the original anndata was also same batch and label
    adata1 = synthetic_iid(run_setup_anndata=False)
    setup_anndata(adata1)
    adata2 = synthetic_iid(run_setup_anndata=False)
    del adata2.obs["batch"]
    transfer_anndata_setup(adata1, adata2)
    assert adata2.obs["_scvi_batch"][0] == 0
    assert adata2.obs["_scvi_labels"][0] == 0


def test_data_format():
    # if data was dense np array, check after setup_anndata, data is C_CONTIGUOUS
    adata = synthetic_iid(run_setup_anndata=False)

    old_x = adata.X
    old_pro = adata.obsm["protein_expression"]
    old_obs = adata.obs
    adata.X = np.asfortranarray(old_x)
    adata.obsm["protein_expression"] = np.asfortranarray(old_pro)
    assert adata.X.flags["C_CONTIGUOUS"] is False
    assert adata.obsm["protein_expression"].flags["C_CONTIGUOUS"] is False

    setup_anndata(adata, protein_expression_obsm_key="protein_expression")
    assert adata.X.flags["C_CONTIGUOUS"] is True
    assert adata.obsm["protein_expression"].flags["C_CONTIGUOUS"] is True

    assert np.array_equal(old_x, adata.X)
    assert np.array_equal(old_pro, adata.obsm["protein_expression"])
    assert np.array_equal(old_obs, adata.obs)

    assert np.array_equal(adata.X, get_from_registry(adata, _CONSTANTS.X_KEY))
    assert np.array_equal(
        adata.obsm["protein_expression"],
        get_from_registry(adata, _CONSTANTS.PROTEIN_EXP_KEY),
    )

    # if data is sparse, check that after setup_anndata, data is csr_matrix
    adata = synthetic_iid(run_setup_anndata=False)
    adata.X = csc_matrix(adata.X)
    adata.obsm["protein_expression"] = csc_matrix(adata.obsm["protein_expression"])
    old_x = adata.X
    old_pro = adata.obsm["protein_expression"]
    setup_anndata(adata, protein_expression_obsm_key="protein_expression")

    assert adata.X.getformat() == "csr"
    assert adata.obsm["protein_expression"].getformat() == "csr"

    assert np.array_equal(old_x.toarray(), adata.X.toarray())
    assert np.array_equal(old_pro.toarray(), adata.obsm["protein_expression"].toarray())

    assert np.array_equal(
        adata.X.toarray(), get_from_registry(adata, _CONSTANTS.X_KEY).toarray()
    )
    assert np.array_equal(
        adata.obsm["protein_expression"].toarray(),
        get_from_registry(adata, _CONSTANTS.PROTEIN_EXP_KEY).toarray(),
    )

    # if obsm is dataframe, make it C_CONTIGUOUS if it isnt
    adata = synthetic_iid()
    pe = np.asfortranarray(adata.obsm["protein_expression"])
    adata.obsm["protein_expression"] = pd.DataFrame(pe, index=adata.obs_names)
    assert adata.obsm["protein_expression"].to_numpy().flags["C_CONTIGUOUS"] is False
    setup_anndata(adata, protein_expression_obsm_key="protein_expression")
    new_pe = get_from_registry(adata, "protein_expression")
    assert new_pe.to_numpy().flags["C_CONTIGUOUS"] is True
    assert np.array_equal(pe, new_pe)
    assert np.array_equal(adata.X, get_from_registry(adata, _CONSTANTS.X_KEY))
    assert np.array_equal(
        adata.obsm["protein_expression"],
        get_from_registry(adata, _CONSTANTS.PROTEIN_EXP_KEY),
    )


def test_setup_anndata():
    # test regular setup
    adata = synthetic_iid(run_setup_anndata=False)
    setup_anndata(
        adata,
        batch_key="batch",
        labels_key="labels",
        protein_expression_obsm_key="protein_expression",
        protein_names_uns_key="protein_names",
    )
    np.testing.assert_array_equal(
        get_from_registry(adata, "batch_indices"),
        np.array(adata.obs["_scvi_batch"]).reshape((-1, 1)),
    )
    np.testing.assert_array_equal(
        get_from_registry(adata, "labels"),
        np.array(adata.obs["labels"].cat.codes).reshape((-1, 1)),
    )
    np.testing.assert_array_equal(get_from_registry(adata, "X"), adata.X)
    np.testing.assert_array_equal(
        get_from_registry(adata, "protein_expression"),
        adata.obsm["protein_expression"],
    )
    np.testing.assert_array_equal(
        adata.uns["scvi_protein_names"], adata.uns["protein_names"]
    )

    # test that error is thrown if its a view:
    adata = synthetic_iid()
    with pytest.raises(ValueError):
        setup_anndata(adata[1])

    # If obsm is a df and protein_names_uns_key is None, protein names should be grabbed from column of df
    adata = synthetic_iid()
    new_protein_names = np.array(random.sample(range(100), 100)).astype("str")
    df = pd.DataFrame(
        adata.obsm["protein_expression"],
        index=adata.obs_names,
        columns=new_protein_names,
    )
    adata.obsm["protein_expression"] = df
    setup_anndata(adata, protein_expression_obsm_key="protein_expression")
    np.testing.assert_array_equal(adata.uns["scvi_protein_names"], new_protein_names)

    # test that layer is working properly
    adata = synthetic_iid()
    true_x = adata.X
    adata.layers["X"] = true_x
    adata.X = np.ones_like(adata.X)
    setup_anndata(adata, layer="X")
    np.testing.assert_array_equal(get_from_registry(adata, "X"), true_x)

    # test that it creates layers and batch if no layers_key is passed
    adata = synthetic_iid()
    setup_anndata(
        adata,
        protein_expression_obsm_key="protein_expression",
        protein_names_uns_key="protein_names",
    )
    np.testing.assert_array_equal(
        get_from_registry(adata, "batch_indices"), np.zeros((adata.shape[0], 1))
    )
    np.testing.assert_array_equal(
        get_from_registry(adata, "labels"), np.zeros((adata.shape[0], 1))
    )


def test_save_setup_anndata(save_path):
    adata = synthetic_iid()
    adata.write(os.path.join(save_path, "test.h5ad"))


def test_extra_covariates():
    adata = synthetic_iid()
    adata.obs["cont1"] = np.random.normal(size=(adata.shape[0],))
    adata.obs["cont2"] = np.random.normal(size=(adata.shape[0],))
    adata.obs["cat1"] = np.random.randint(0, 5, size=(adata.shape[0],))
    adata.obs["cat2"] = np.random.randint(0, 5, size=(adata.shape[0],))
    setup_anndata(
        adata,
        batch_key="batch",
        labels_key="labels",
        protein_expression_obsm_key="protein_expression",
        protein_names_uns_key="protein_names",
        continuous_covariate_keys=["cont1", "cont2"],
        categorical_covariate_keys=["cat1", "cat2"],
    )
    df1 = adata.obsm["_scvi_extra_continuous"]
    df2 = adata.obs[["cont1", "cont2"]]
    pd.testing.assert_frame_equal(df1, df2)


def test_register_tensor_from_anndata():
    adata = synthetic_iid()
    adata.obs["cont1"] = np.random.normal(size=(adata.shape[0],))
    register_tensor_from_anndata(
        adata, registry_key="test", adata_attr_name="obs", adata_key_name="cont1"
    )
    assert "test" in adata.uns["_scvi"]["data_registry"]
    assert adata.uns["_scvi"]["data_registry"]["test"] == dict(
        attr_name="obs", attr_key="cont1"
    )


def test_scvidataset_getitem():
    adata = synthetic_iid()
    setup_anndata(
        adata,
        batch_key="batch",
        labels_key="labels",
        protein_expression_obsm_key="protein_expression",
        protein_names_uns_key="protein_names",
    )
    # check that we can successfully pass in a list of tensors to get
    tensors_to_get = ["batch_indices", "local_l_var"]
    bd = ScviDataset(adata, getitem_tensors=tensors_to_get)
    np.testing.assert_array_equal(tensors_to_get, list(bd[1].keys()))

    # check that we can successfully pass in a dict of tensors and their associated types
    bd = ScviDataset(adata, getitem_tensors={"X": np.int, "local_l_var": np.float64})
    assert bd[1]["X"].dtype == np.int64
    assert bd[1]["local_l_var"].dtype == np.float64

    # check that by default we get all the registered tensors
    bd = ScviDataset(adata)
    all_registered_tensors = list(adata.uns["_scvi"]["data_registry"].keys())
    np.testing.assert_array_equal(all_registered_tensors, list(bd[1].keys()))
    assert bd[1]["X"].shape[0] == bd.n_vars

    # check that ScviDataset returns numpy array
    adata1 = synthetic_iid()
    bd = ScviDataset(adata1)
    for key, value in bd[1].items():
        assert type(value) == np.ndarray

    # check ScviDataset returns numpy array counts were sparse
    adata = synthetic_iid(run_setup_anndata=False)
    adata.X = sparse.csr_matrix(adata.X)
    setup_anndata(adata)
    bd = ScviDataset(adata)
    for key, value in bd[1].items():
        assert type(value) == np.ndarray

    # check ScviDataset returns numpy array if pro exp was sparse
    adata = synthetic_iid(run_setup_anndata=False)
    adata.obsm["protein_expression"] = sparse.csr_matrix(
        adata.obsm["protein_expression"]
    )
    setup_anndata(
        adata, batch_key="batch", protein_expression_obsm_key="protein_expression"
    )
    bd = ScviDataset(adata)
    for key, value in bd[1].items():
        assert type(value) == np.ndarray

    # check pro exp is being returned as numpy array even if its DF
    adata = synthetic_iid(run_setup_anndata=False)
    adata.obsm["protein_expression"] = pd.DataFrame(
        adata.obsm["protein_expression"], index=adata.obs_names
    )
    setup_anndata(
        adata, batch_key="batch", protein_expression_obsm_key="protein_expression"
    )
    bd = ScviDataset(adata)
    for key, value in bd[1].items():
        assert type(value) == np.ndarray
