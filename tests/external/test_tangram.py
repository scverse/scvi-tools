import mudata

from scvi.data import synthetic_iid
from scvi.external import Tangram


def _get_mdata():
    dataset1 = synthetic_iid(n_labels=5, batch_size=100)
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
    model.train(max_epochs=1)
