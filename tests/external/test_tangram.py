import mudata

from scvi.data import synthetic_iid
from scvi.external import Tangram


def _get_mdata(sparse=False):
    dataset1 = synthetic_iid(batch_size=100, sparse=sparse)
    dataset2 = dataset1[-25:].copy()
    dataset1 = dataset1[:-25].copy()
    mdata = mudata.MuData({"sc": dataset1, "sp": dataset2})
    return mdata


def test_tangram():
    mdata = _get_mdata()
    Tangram.setup_mudata(
        mdata,
        density_prior_key=None,
        modalities={"density_prior_key": "sp", "sc_layer": "sc", "sp_layer": "sp"},
    )
    model = Tangram(mdata)
    model.train(max_epochs=1, retain_sparsity=False)


def test_tangram_sparse_input():
    mdata = _get_mdata(sparse=True)
    Tangram.setup_mudata(
        mdata,
        density_prior_key=None,
        modalities={"density_prior_key": "sp", "sc_layer": "sc", "sp_layer": "sp"},
    )
    model = Tangram(mdata)
    model.train(max_epochs=1)
