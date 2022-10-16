import mudata
import numpy as np
import pytest

from scvi.data import synthetic_iid
from scvi.external import Tangram


def _get_mdata(sparse=False):
    dataset1 = synthetic_iid(batch_size=100, sparse=sparse)
    dataset2 = dataset1[-25:].copy()
    dataset1 = dataset1[:-25].copy()
    mdata = mudata.MuData({"sc": dataset1, "sp": dataset2})
    ad_sp = mdata.mod["sp"]
    rna_count_per_spot = np.asarray(ad_sp.X.sum(axis=1)).squeeze()
    ad_sp.obs["rna_count_based_density"] = rna_count_per_spot / np.sum(
        rna_count_per_spot
    )
    return mdata


@pytest.mark.parametrize("density_prior_key", [None, "rna_count_based_density"])
def test_tangram(density_prior_key):
    mdata = _get_mdata()
    Tangram.setup_mudata(
        mdata,
        density_prior_key=density_prior_key,
        modalities={"density_prior_key": "sp", "sc_layer": "sc", "sp_layer": "sp"},
    )
    model = Tangram(mdata)
    model.train(max_epochs=1)
    model.get_mapper_matrix()


def test_tangram_constrained():
    mdata = _get_mdata()
    Tangram.setup_mudata(
        mdata,
        density_prior_key="rna_count_based_density",
        modalities={"density_prior_key": "sp", "sc_layer": "sc", "sp_layer": "sp"},
    )
    model = Tangram(
        mdata,
        constrained=True,
        target_count=2,
    )
    model.train(max_epochs=1)
    model.get_mapper_matrix()
