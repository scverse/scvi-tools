import mudata
import numpy as np

from scvi import REGISTRY_KEYS
from scvi.data import synthetic_iid

from .utils import generic_setup_mudata_manager


def test_setup_mudata():
    # test regular setup
    adata = synthetic_iid()
    mdata = mudata.MuData({"rna": adata})
    adata_manager = generic_setup_mudata_manager(mdata, layer=None, layer_mod="rna")
    np.testing.assert_array_equal(
        adata_manager.get_from_registry(REGISTRY_KEYS.X_KEY), adata.X
    )
