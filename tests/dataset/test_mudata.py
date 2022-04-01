import mudata
import numpy as np
import pytest

from scvi import REGISTRY_KEYS
from scvi.data import synthetic_iid

from .utils import generic_setup_mudata_manager


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
