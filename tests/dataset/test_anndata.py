import os
import random
from typing import List, Optional

import anndata
import numpy as np
import pandas as pd
import pytest
import scipy.sparse as sparse
from scipy.sparse.csr import csr_matrix

import scvi
from scvi import _CONSTANTS
from scvi.data import synthetic_iid
from scvi.data.anndata import AnnDataManager
from scvi.data.anndata.fields import (
    CategoricalJointObsField,
    CategoricalObsField,
    LayerField,
    NumericalJointObsField,
    ProteinObsmField,
)
from scvi.dataloaders import AnnTorchDataset


def _generic_setup_adata_manager(
    adata: anndata.AnnData,
    batch_key: Optional[str] = None,
    labels_key: Optional[str] = None,
    categorical_covariate_keys: Optional[List[str]] = None,
    continuous_covariate_keys: Optional[List[str]] = None,
    layer: Optional[str] = None,
    protein_expression_obsm_key: Optional[str] = None,
    protein_names_uns_key: Optional[str] = None,
) -> AnnDataManager:
    batch_field = CategoricalObsField(_CONSTANTS.BATCH_KEY, batch_key)
    anndata_fields = [
        batch_field,
        LayerField(_CONSTANTS.X_KEY, layer, is_count_data=True),
        CategoricalObsField(_CONSTANTS.LABELS_KEY, labels_key),
        CategoricalJointObsField(_CONSTANTS.CAT_COVS_KEY, categorical_covariate_keys),
        NumericalJointObsField(_CONSTANTS.CONT_COVS_KEY, continuous_covariate_keys),
    ]
    if protein_expression_obsm_key is not None:
        anndata_fields.append(
            ProteinObsmField(
                _CONSTANTS.PROTEIN_EXP_KEY,
                protein_expression_obsm_key,
                batch_field.attr_key,
                colnames_uns_key=protein_names_uns_key,
                is_count_data=True,
            )
        )
    adata_manager = AnnDataManager(fields=anndata_fields)
    adata_manager.register_fields(adata)
    return adata_manager


def test_transfer_anndata_setup():
    # test transfer_anndata function
    adata1 = synthetic_iid()
    adata2 = synthetic_iid()
    adata2.X = adata1.X
    adata1_manager = _generic_setup_adata_manager(adata1)
    adata1_manager.transfer_setup(adata2)
    np.testing.assert_array_equal(
        adata1.obs["_scvi_labels"], adata2.obs["_scvi_labels"]
    )

    # test if layer was used initially, again used in transfer setup
    adata1 = synthetic_iid()
    adata2 = synthetic_iid()
    raw_counts = adata1.X.copy()
    adata1.layers["raw"] = raw_counts
    adata2.layers["raw"] = raw_counts
    zeros = np.zeros_like(adata1.X)
    ones = np.ones_like(adata1.X)
    adata1.X = zeros
    adata2.X = ones
    adata1_manager = _generic_setup_adata_manager(adata1, layer="raw")
    adata1_manager.transfer_setup(adata2)
    np.testing.assert_array_equal(
        adata1.obs["_scvi_labels"], adata2.obs["_scvi_labels"]
    )

    # test that an unknown batch throws an error
    adata1 = synthetic_iid()
    adata2 = synthetic_iid()
    adata2.obs["batch"] = [2] * adata2.n_obs
    adata1_manager = _generic_setup_adata_manager(adata1, batch_key="batch")
    with pytest.raises(ValueError):
        adata1_manager.transfer_setup(adata2)

    # TODO: test that a batch with wrong dtype throws an error
    # adata1 = synthetic_iid()
    # adata2 = synthetic_iid()
    # adata2.obs["batch"] = ["0"] * adata2.n_obs
    # with pytest.raises(ValueError):
    #     transfer_anndata_setup(adata1, adata2)

    # test that an unknown label throws an error
    adata1 = synthetic_iid()
    adata2 = synthetic_iid()
    adata2.obs["labels"] = ["label_123"] * adata2.n_obs
    adata1_manager = _generic_setup_adata_manager(adata1, labels_key="labels")
    with pytest.raises(ValueError):
        adata1_manager.transfer_setup(adata2)

    # test that correct mapping was applied
    adata1 = synthetic_iid()
    adata2 = synthetic_iid()
    adata2.obs["labels"] = ["label_1"] * adata2.n_obs
    adata1_manager = _generic_setup_adata_manager(adata1, labels_key="labels")
    adata1_manager.transfer_setup(adata2)
    labels_mapping = adata1_manager.get_state_registry("labels")[
        CategoricalObsField.CATEGORICAL_MAPPING_KEY
    ]
    correct_label = np.where(labels_mapping == "label_1")[0][0]
    adata2.obs["_scvi_labels"][0] == correct_label

    # test that transfer_anndata_setup correctly looks for adata.obs['batch']
    adata1 = synthetic_iid()
    adata2 = synthetic_iid()
    del adata2.obs["batch"]
    adata1_manager = _generic_setup_adata_manager(adata1, batch_key="batch")
    with pytest.raises(KeyError):
        adata1_manager.transfer_setup(adata2)

    # test that transfer_anndata_setup assigns same batch and label to cells
    # if the original anndata was also same batch and label
    adata1 = synthetic_iid()
    adata1_manager = _generic_setup_adata_manager(adata1)
    adata2 = synthetic_iid()
    del adata2.obs["batch"]
    adata1_manager.transfer_setup(adata2)
    assert adata2.obs["_scvi_batch"][0] == 0
    assert adata2.obs["_scvi_labels"][0] == 0

    # test that if a category mapping is a subset, transfer anndata is called
    a1 = scvi.data.synthetic_iid()
    scvi.model.SCVI.setup_anndata(a1, batch_key="batch")
    a2 = scvi.data.synthetic_iid()
    a2.obs["batch"] = "batch_1"
    scvi.model.SCVI.setup_anndata(a2, batch_key="batch")
    m = scvi.model.SCVI(a1)
    m.train(1)
    m.get_latent_representation(a2)
    assert a2.obs["_scvi_batch"].all() == 1


def test_data_format():
    # if data was dense np array, check after setup_anndata, data is C_CONTIGUOUS
    adata = synthetic_iid()

    old_x = adata.X
    old_pro = adata.obsm["protein_expression"]
    old_obs = adata.obs
    adata.X = np.asfortranarray(old_x)
    adata.obsm["protein_expression"] = np.asfortranarray(old_pro)
    assert adata.X.flags["C_CONTIGUOUS"] is False
    assert adata.obsm["protein_expression"].flags["C_CONTIGUOUS"] is False

    adata_manager = _generic_setup_adata_manager(
        adata, protein_expression_obsm_key="protein_expression"
    )
    assert adata.X.flags["C_CONTIGUOUS"] is True
    assert adata.obsm["protein_expression"].flags["C_CONTIGUOUS"] is True

    assert np.array_equal(old_x, adata.X)
    assert np.array_equal(old_pro, adata.obsm["protein_expression"])
    assert np.array_equal(old_obs, adata.obs)

    assert np.array_equal(adata.X, adata_manager.get_from_registry(_CONSTANTS.X_KEY))
    assert np.array_equal(
        adata.obsm["protein_expression"],
        adata_manager.get_from_registry(_CONSTANTS.PROTEIN_EXP_KEY),
    )

    # if obsm is dataframe, make it C_CONTIGUOUS if it isnt
    adata = synthetic_iid()
    pe = np.asfortranarray(adata.obsm["protein_expression"])
    adata.obsm["protein_expression"] = pd.DataFrame(pe, index=adata.obs_names)
    assert adata.obsm["protein_expression"].to_numpy().flags["C_CONTIGUOUS"] is False
    adata_manager = _generic_setup_adata_manager(
        adata, protein_expression_obsm_key="protein_expression"
    )
    new_pe = adata_manager.get_from_registry(_CONSTANTS.PROTEIN_EXP_KEY)
    assert new_pe.to_numpy().flags["C_CONTIGUOUS"] is True
    assert np.array_equal(pe, new_pe)
    assert np.array_equal(adata.X, adata_manager.get_from_registry(_CONSTANTS.X_KEY))
    assert np.array_equal(
        adata.obsm["protein_expression"],
        adata_manager.get_from_registry(_CONSTANTS.PROTEIN_EXP_KEY),
    )


def test_setup_anndata():
    # test regular setup
    adata = synthetic_iid()
    adata_manager = _generic_setup_adata_manager(
        adata,
        batch_key="batch",
        labels_key="labels",
        protein_expression_obsm_key="protein_expression",
        protein_names_uns_key="protein_names",
    )
    np.testing.assert_array_equal(
        adata_manager.get_from_registry(_CONSTANTS.BATCH_KEY),
        np.array(adata.obs["_scvi_batch"]).reshape((-1, 1)),
    )
    np.testing.assert_array_equal(
        adata_manager.get_from_registry(_CONSTANTS.LABELS_KEY),
        np.array(adata.obs["labels"].cat.codes).reshape((-1, 1)),
    )
    np.testing.assert_array_equal(
        adata_manager.get_from_registry(_CONSTANTS.X_KEY), adata.X
    )
    np.testing.assert_array_equal(
        adata_manager.get_from_registry(_CONSTANTS.PROTEIN_EXP_KEY),
        adata.obsm["protein_expression"],
    )
    np.testing.assert_array_equal(
        adata_manager.get_state_registry(_CONSTANTS.PROTEIN_EXP_KEY)[
            ProteinObsmField.COLUMN_NAMES_KEY
        ],
        adata.uns["protein_names"],
    )

    # test that error is thrown if its a view:
    adata = synthetic_iid()
    with pytest.raises(ValueError):
        _generic_setup_adata_manager(adata[1])

    # If obsm is a df and protein_names_uns_key is None, protein names should be grabbed from column of df
    adata = synthetic_iid()
    new_protein_names = np.array(random.sample(range(100), 100)).astype("str")
    df = pd.DataFrame(
        adata.obsm["protein_expression"],
        index=adata.obs_names,
        columns=new_protein_names,
    )
    adata.obsm["protein_expression"] = df
    adata_manager = _generic_setup_adata_manager(
        adata, protein_expression_obsm_key="protein_expression"
    )
    np.testing.assert_array_equal(
        adata_manager.get_state_registry(_CONSTANTS.PROTEIN_EXP_KEY)[
            ProteinObsmField.COLUMN_NAMES_KEY
        ],
        new_protein_names,
    )

    # test that layer is working properly
    adata = synthetic_iid()
    true_x = adata.X
    adata.layers["X"] = true_x
    adata.X = np.ones_like(adata.X)
    adata_manager = _generic_setup_adata_manager(adata, layer="X")
    np.testing.assert_array_equal(
        adata_manager.get_from_registry(_CONSTANTS.X_KEY), true_x
    )

    # test that it creates layers and batch if no layers_key is passed
    adata = synthetic_iid()
    adata_manager = _generic_setup_adata_manager(
        adata,
        protein_expression_obsm_key="protein_expression",
        protein_names_uns_key="protein_names",
    )
    np.testing.assert_array_equal(
        adata_manager.get_from_registry(_CONSTANTS.BATCH_KEY),
        np.zeros((adata.shape[0], 1)),
    )
    np.testing.assert_array_equal(
        adata_manager.get_from_registry(_CONSTANTS.LABELS_KEY),
        np.zeros((adata.shape[0], 1)),
    )


def test_save_setup_anndata(save_path):
    adata = synthetic_iid()
    _generic_setup_adata_manager(
        adata,
        batch_key="batch",
        labels_key="labels",
        protein_expression_obsm_key="protein_expression",
        protein_names_uns_key="protein_names",
    )
    adata.write(os.path.join(save_path, "test.h5ad"))


def test_extra_covariates():
    adata = synthetic_iid()
    adata.obs["cont1"] = np.random.normal(size=(adata.shape[0],))
    adata.obs["cont2"] = np.random.normal(size=(adata.shape[0],))
    adata.obs["cat1"] = np.random.randint(0, 5, size=(adata.shape[0],))
    adata.obs["cat2"] = np.random.randint(0, 5, size=(adata.shape[0],))
    scvi.model.SCVI.setup_anndata(
        adata,
        batch_key="batch",
        labels_key="labels",
        continuous_covariate_keys=["cont1", "cont2"],
        categorical_covariate_keys=["cat1", "cat2"],
    )
    m = scvi.model.SCVI(adata)
    m.train(1)
    df1 = m.get_from_registry(adata, _CONSTANTS.CONT_COVS_KEY)
    df2 = adata.obs[["cont1", "cont2"]]
    pd.testing.assert_frame_equal(df1, df2)


def test_extra_covariates_transfer():
    adata = synthetic_iid()
    adata.obs["cont1"] = np.random.normal(size=(adata.shape[0],))
    adata.obs["cont2"] = np.random.normal(size=(adata.shape[0],))
    adata.obs["cat1"] = np.random.randint(0, 5, size=(adata.shape[0],))
    adata.obs["cat2"] = np.random.randint(0, 5, size=(adata.shape[0],))
    adata_manager = _generic_setup_adata_manager(
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

    adata_manager.transfer_setup(bdata)

    # give it a new category
    bdata.obs["cat1"] = 6
    bdata_manager = adata_manager.transfer_setup(bdata, extend_categories=True)
    assert (
        bdata_manager.get_state_registry(_CONSTANTS.CAT_COVS_KEY)[
            CategoricalJointObsField.MAPPINGS_KEY
        ]["cat1"][-1]
        == 6
    )


def test_anntorchdataset_getitem():
    adata = synthetic_iid()
    adata_manager = _generic_setup_adata_manager(
        adata,
        batch_key="batch",
        labels_key="labels",
        protein_expression_obsm_key="protein_expression",
        protein_names_uns_key="protein_names",
    )
    # check that we can successfully pass in a list of tensors to get
    tensors_to_get = [_CONSTANTS.BATCH_KEY, _CONSTANTS.LABELS_KEY]
    bd = AnnTorchDataset(adata_manager, getitem_tensors=tensors_to_get)
    np.testing.assert_array_equal(tensors_to_get, list(bd[1].keys()))

    # check that we can successfully pass in a dict of tensors and their associated types
    bd = AnnTorchDataset(
        adata_manager,
        getitem_tensors={_CONSTANTS.X_KEY: np.int, _CONSTANTS.LABELS_KEY: np.int64},
    )
    assert bd[1][_CONSTANTS.X_KEY].dtype == np.int64
    assert bd[1][_CONSTANTS.LABELS_KEY].dtype == np.int64

    # check that by default we get all the registered tensors
    bd = AnnTorchDataset(adata_manager)
    all_registered_tensors = list(adata_manager.data_registry.keys())
    np.testing.assert_array_equal(all_registered_tensors, list(bd[1].keys()))
    assert bd[1][_CONSTANTS.X_KEY].shape[0] == bd.adata_manager.summary_stats["n_vars"]

    # check that AnnTorchDataset returns numpy array
    adata1 = synthetic_iid()
    adata1_manager = _generic_setup_adata_manager(adata1)
    bd = AnnTorchDataset(adata1_manager)
    for value in bd[1].values():
        assert type(value) == np.ndarray

    # check AnnTorchDataset returns numpy array counts were sparse
    adata = synthetic_iid()
    adata.X = sparse.csr_matrix(adata.X)
    adata_manager = _generic_setup_adata_manager(adata1)
    bd = AnnTorchDataset(adata_manager)
    for value in bd[1].values():
        assert type(value) == np.ndarray

    # check AnnTorchDataset returns numpy array if pro exp was sparse
    adata = synthetic_iid()
    adata.obsm["protein_expression"] = sparse.csr_matrix(
        adata.obsm["protein_expression"]
    )
    adata_manager = _generic_setup_adata_manager(
        adata, batch_key="batch", protein_expression_obsm_key="protein_expression"
    )
    bd = AnnTorchDataset(adata_manager)
    for value in bd[1].values():
        assert type(value) == np.ndarray

    # check pro exp is being returned as numpy array even if its DF
    adata = synthetic_iid()
    adata.obsm["protein_expression"] = pd.DataFrame(
        adata.obsm["protein_expression"], index=adata.obs_names
    )
    _generic_setup_adata_manager(
        adata, batch_key="batch", protein_expression_obsm_key="protein_expression"
    )
    bd = AnnTorchDataset(adata_manager)
    for value in bd[1].values():
        assert type(value) == np.ndarray


def test_saving(save_path):
    save_path = os.path.join(save_path, "tmp_adata.h5ad")
    adata = synthetic_iid()
    adata.obs["cont1"] = np.random.uniform(5, adata.n_obs)
    adata.obs["cont2"] = np.random.uniform(5, adata.n_obs)
    adata.obs["cat1"] = np.random.randint(0, 3, adata.n_obs).astype(str)
    adata.obs["cat1"][1] = "asdf"
    adata.obs["cat1"][2] = "f34"
    adata.obs["cat2"] = np.random.randint(0, 7, adata.n_obs)

    _generic_setup_adata_manager(
        adata,
        protein_expression_obsm_key="protein_expression",
        batch_key="batch",
        labels_key="labels",
        categorical_covariate_keys=["cat1", "cat2"],
        continuous_covariate_keys=["cont1", "cont2"],
    )
    adata.write(save_path)
    anndata.read(save_path)


def test_backed_anndata(save_path):
    adata = scvi.data.synthetic_iid()
    path = os.path.join(save_path, "test_data.h5ad")
    adata.write_h5ad(path)
    adata = anndata.read_h5ad(path, backed="r+")
    adata_manager = _generic_setup_adata_manager(adata, batch_key="batch")

    # test get item
    bd = AnnTorchDataset(adata_manager)
    bd[np.arange(adata.n_obs)]

    # sparse
    adata = scvi.data.synthetic_iid()
    adata.X = csr_matrix(adata.X)
    path = os.path.join(save_path, "test_data2.h5ad")
    adata.write_h5ad(path)
    adata = anndata.read_h5ad(path, backed="r+")
    adata_manager = _generic_setup_adata_manager(adata, batch_key="batch")

    # test get item
    bd = AnnTorchDataset(adata_manager)
    bd[np.arange(adata.n_obs)]
