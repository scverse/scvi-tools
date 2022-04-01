import os

import mudata
import numpy as np
import pytest

from scvi import REGISTRY_KEYS
from scvi.data import synthetic_iid

from .utils import generic_setup_mudata_manager


def test_transfer_fields():
    # test transfer_fields function
    adata1 = synthetic_iid()
    protein_adata1 = synthetic_iid()
    mdata1 = mudata.MuData({"rna": adata1, "protein": protein_adata1})
    adata2 = synthetic_iid()
    adata2.X = adata1.X
    protein_adata2 = synthetic_iid()
    mdata2 = mudata.MuData({"rna": adata2, "protein": protein_adata2})

    adata1_manager = generic_setup_mudata_manager(
        mdata1,
        layer_mod="rna",
        layer=None,
        batch_mod="rna",
        batch_key="batch",
        protein_expression_mod="protein",
        protein_expression_layer=None,
    )
    adata1_manager.transfer_fields(mdata2)
    np.testing.assert_array_equal(adata1.obs["_scvi_batch"], adata2.obs["_scvi_batch"])


def test_transfer_fields_view():
    # test transfer_fields function
    adata1 = synthetic_iid()
    protein_adata1 = synthetic_iid()
    mdata1 = mudata.MuData({"rna": adata1, "protein": protein_adata1})
    adata2 = synthetic_iid()
    adata2.X = adata1.X
    protein_adata2 = synthetic_iid()
    mdata2 = mudata.MuData({"rna": adata2, "protein": protein_adata2})[:50]

    adata1_manager = generic_setup_mudata_manager(
        mdata1,
        layer_mod="rna",
        layer=None,
        batch_mod="rna",
        batch_key="batch",
        protein_expression_mod="protein",
        protein_expression_layer=None,
    )
    with pytest.raises(ValueError):
        adata1_manager.transfer_fields(mdata2)


def test_transfer_fields_layer():
    # test if layer was used initially, again used in transfer setup
    adata1 = synthetic_iid()
    raw_counts = adata1.X.copy()
    adata1.layers["raw"] = raw_counts
    protein_adata1 = synthetic_iid()
    mdata1 = mudata.MuData({"rna": adata1, "protein": protein_adata1})
    adata2 = synthetic_iid()
    adata2.layers["raw"] = raw_counts
    protein_adata2 = synthetic_iid()
    mdata2 = mudata.MuData({"rna": adata2, "protein": protein_adata2})

    zeros = np.zeros_like(adata1.X)
    ones = np.ones_like(adata1.X)
    adata1.X = zeros
    adata2.X = ones

    adata_manager = generic_setup_mudata_manager(
        mdata1,
        layer_mod="rna",
        layer="raw",
        batch_mod="rna",
        batch_key="batch",
        protein_expression_mod="protein",
        protein_expression_layer=None,
    )
    adata_manager.transfer_fields(mdata2)
    np.testing.assert_array_equal(adata1.obs["_scvi_batch"], adata2.obs["_scvi_batch"])


def test_transfer_fields_unknown_batch():
    # test that an unknown batch throws an error
    adata1 = synthetic_iid()
    protein_adata1 = synthetic_iid()
    mdata1 = mudata.MuData({"rna": adata1, "protein": protein_adata1})
    adata2 = synthetic_iid()
    adata2.X = adata1.X
    protein_adata2 = synthetic_iid()
    mdata2 = mudata.MuData({"rna": adata2, "protein": protein_adata2})
    adata2.obs["batch"] = [2] * adata2.n_obs
    adata1_manager = generic_setup_mudata_manager(
        mdata1, layer_mod="rna", batch_mod="rna", batch_key="batch"
    )
    with pytest.raises(ValueError):
        adata1_manager.transfer_fields(mdata2)


def test_transfer_fields_diff_batch_mapping():
    # test that correct mapping was applied
    adata1 = synthetic_iid()
    protein_adata1 = synthetic_iid()
    mdata1 = mudata.MuData({"rna": adata1, "protein": protein_adata1})
    adata2 = synthetic_iid()
    adata2.X = adata1.X
    protein_adata2 = synthetic_iid()
    mdata2 = mudata.MuData({"rna": adata2, "protein": protein_adata2})
    adata2.obs["labels"] = ["batch_1"] * adata2.n_obs
    adata1_manager = generic_setup_mudata_manager(
        mdata1, layer_mod="rna", batch_mod="rna", batch_key="batch"
    )
    adata1_manager.transfer_fields(mdata2)
    batch_mapping = adata1_manager.get_state_registry(
        REGISTRY_KEYS.BATCH_KEY
    ).categorical_mapping
    correct_batch = np.where(batch_mapping == "batch_1")[0][0]
    adata2.obs["_scvi_batch"][0] == correct_batch


def test_transfer_fields_missing_batch():
    # test that transfer_fields correctly looks for adata.obs['batch']
    adata1 = synthetic_iid()
    protein_adata1 = synthetic_iid()
    mdata1 = mudata.MuData({"rna": adata1, "protein": protein_adata1})
    adata2 = synthetic_iid()
    del adata2.obs["batch"]
    adata2.X = adata1.X
    protein_adata2 = synthetic_iid()
    mdata2 = mudata.MuData({"rna": adata2, "protein": protein_adata2})
    adata1_manager = generic_setup_mudata_manager(
        mdata1, layer_mod="rna", batch_mod="rna", batch_key="batch"
    )
    with pytest.raises(KeyError):
        adata1_manager.transfer_fields(mdata2)


def test_transfer_fields_default_batch():
    # test that transfer_fields assigns same batch to cells
    # if the original anndata was also same batch
    adata1 = synthetic_iid()
    protein_adata1 = synthetic_iid()
    mdata1 = mudata.MuData({"rna": adata1, "protein": protein_adata1})
    adata2 = synthetic_iid()
    del adata2.obs["batch"]
    adata2.X = adata1.X
    protein_adata2 = synthetic_iid()
    mdata2 = mudata.MuData({"rna": adata2, "protein": protein_adata2})
    adata1_manager = generic_setup_mudata_manager(
        mdata1, layer_mod="rna", batch_mod="rna", batch_key=None
    )
    adata1_manager.transfer_fields(mdata2)
    assert adata2.obs["_scvi_batch"][0] == 0


def test_data_format():
    # if data was dense np array, check after setup_anndata, data is C_CONTIGUOUS
    adata = synthetic_iid()
    protein_adata = synthetic_iid()
    mdata = mudata.MuData({"rna": adata, "protein": protein_adata})

    old_x = adata.X
    old_pro = protein_adata.X
    old_obs = adata.obs
    adata.X = np.asfortranarray(old_x)
    protein_adata.X = np.asfortranarray(old_pro)
    assert adata.X.flags["C_CONTIGUOUS"] is False
    assert protein_adata.X.flags["C_CONTIGUOUS"] is False

    adata_manager = generic_setup_mudata_manager(
        mdata, layer_mod="rna", protein_expression_mod="protein"
    )
    assert adata.X.flags["C_CONTIGUOUS"] is True
    assert protein_adata.X.flags["C_CONTIGUOUS"] is True

    assert np.array_equal(old_x, adata.X)
    assert np.array_equal(old_pro, protein_adata.X)
    assert np.array_equal(old_obs, adata.obs)

    assert np.array_equal(adata.X, adata_manager.get_from_registry(REGISTRY_KEYS.X_KEY))
    assert np.array_equal(
        protein_adata.X,
        adata_manager.get_from_registry(REGISTRY_KEYS.PROTEIN_EXP_KEY),
    )

    # if obsm is dataframe, make it C_CONTIGUOUS if it isnt
    adata = synthetic_iid()
    protein_adata = synthetic_iid()
    mdata = mudata.MuData({"rna": adata, "protein": protein_adata})

    pe = np.asfortranarray(protein_adata.X)
    protein_adata.X = pe
    assert protein_adata.X.flags["C_CONTIGUOUS"] is False

    adata_manager = generic_setup_mudata_manager(
        mdata, layer_mod="rna", protein_expression_mod="protein"
    )
    new_pe = adata_manager.get_from_registry(REGISTRY_KEYS.PROTEIN_EXP_KEY)
    assert new_pe.flags["C_CONTIGUOUS"] is True
    assert np.array_equal(pe, new_pe)
    assert np.array_equal(adata.X, adata_manager.get_from_registry(REGISTRY_KEYS.X_KEY))
    assert np.array_equal(
        protein_adata.X,
        adata_manager.get_from_registry(REGISTRY_KEYS.PROTEIN_EXP_KEY),
    )
    assert adata.X.flags["C_CONTIGUOUS"] is True
    assert protein_adata.X.flags["C_CONTIGUOUS"] is True


def test_setup_mudata():
    adata = synthetic_iid()
    protein_adata = synthetic_iid()
    mdata = mudata.MuData({"rna": adata, "protein": protein_adata})
    adata_manager = generic_setup_mudata_manager(
        mdata,
        layer_mod="rna",
        layer=None,
        batch_mod="rna",
        batch_key="batch",
        protein_expression_mod="protein",
        protein_expression_layer=None,
    )
    np.testing.assert_array_equal(
        adata_manager.get_from_registry(REGISTRY_KEYS.X_KEY), adata.X
    )
    np.testing.assert_array_equal(
        adata_manager.get_from_registry(REGISTRY_KEYS.BATCH_KEY),
        np.array(adata.obs["_scvi_batch"]).reshape((-1, 1)),
    )
    np.testing.assert_array_equal(
        adata_manager.get_from_registry(REGISTRY_KEYS.PROTEIN_EXP_KEY),
        protein_adata.X,
    )

    # test that error is thrown if its a view:
    adata = synthetic_iid()
    mdata = mudata.MuData({"rna": adata[1]})
    with pytest.raises(ValueError):
        generic_setup_mudata_manager(mdata, layer_mod="rna")

    adata = synthetic_iid()
    mdata = mudata.MuData({"rna": adata})
    with pytest.raises(ValueError):
        generic_setup_mudata_manager(mdata[1], layer_mod="rna")

    # test that error is thrown if an anndata is passed in:
    adata = synthetic_iid()
    with pytest.raises(AssertionError):
        generic_setup_mudata_manager(adata, layer_mod="rna")

    # test that layer is working properly
    adata = synthetic_iid()
    true_x = adata.X
    adata.layers["X"] = true_x
    adata.X = np.ones_like(adata.X)
    mdata = mudata.MuData({"rna": adata})
    adata_manager = generic_setup_mudata_manager(mdata, layer_mod="rna", layer="X")
    np.testing.assert_array_equal(
        adata_manager.get_from_registry(REGISTRY_KEYS.X_KEY), true_x
    )

    # test that it creates batch if no layers_key is passed
    adata = synthetic_iid()
    mdata = mudata.MuData({"rna": adata})
    adata_manager = generic_setup_mudata_manager(
        mdata, layer_mod="rna", batch_mod="rna"
    )
    np.testing.assert_array_equal(
        adata_manager.get_from_registry(REGISTRY_KEYS.BATCH_KEY),
        np.zeros((adata.shape[0], 1)),
    )
    np.testing.assert_array_equal(
        adata_manager.get_from_registry(REGISTRY_KEYS.BATCH_KEY),
        adata.obs["_scvi_batch"].to_numpy().reshape(-1, 1),
    )

    # test error is thrown when categorical obs field contains nans
    adata = synthetic_iid()
    adata.obs["batch"][:10] = np.nan
    mdata = mudata.MuData({"rna": adata})
    with pytest.raises(ValueError):
        generic_setup_mudata_manager(
            mdata, layer_mod="rna", batch_mod="rna", batch_key="batch"
        )


def test_save_setup_mudata(save_path):
    adata = synthetic_iid()
    mdata = mudata.MuData({"rna": adata})
    generic_setup_mudata_manager(
        mdata,
        layer_mod="rna",
        batch_mod="rna",
        batch_key="batch",
    )
    temp_path = os.path.join(save_path, "test.h5mu")
    mdata.write(temp_path)
    mudata.read(temp_path)


def test_view_registry():
    adata = synthetic_iid()
    protein_adata = synthetic_iid()
    mdata = mudata.MuData({"rna": adata, "protein": protein_adata})
    adata_manager = generic_setup_mudata_manager(
        mdata,
        layer_mod="rna",
        layer=None,
        batch_mod="rna",
        batch_key="batch",
        protein_expression_mod="protein",
        protein_expression_layer=None,
    )
    adata_manager.view_registry()
    adata_manager.view_registry(hide_state_registries=True)
