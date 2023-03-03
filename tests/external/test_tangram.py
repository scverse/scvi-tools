import mudata
import numpy as np
import pytest

from scvi.data import synthetic_iid
from scvi.external import Tangram

modalities = {"density_prior_key": "sp", "sc_layer": "sc", "sp_layer": "sp"}


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
    ad_sp.obs["bad_prior"] = np.random.uniform(size=ad_sp.n_obs)
    return mdata


@pytest.mark.parametrize(
    "density_prior_key,constrained",
    [
        (None, False),
        ("rna_count_based_density", False),
        ("rna_count_based_density", True),
    ],
)
def test_tangram(density_prior_key, constrained):
    mdata = _get_mdata()
    Tangram.setup_mudata(
        mdata,
        density_prior_key=density_prior_key,
        modalities=modalities,
    )
    if constrained:
        target_count = 2
    else:
        target_count = None
    model = Tangram(mdata, constrained=constrained, target_count=target_count)
    model.train(max_epochs=1)
    mdata.mod["sc"].obsm["mapper"] = model.get_mapper_matrix()
    model.project_cell_annotations(
        mdata.mod["sc"],
        mdata.mod["sp"],
        mdata.mod["sc"].obsm["mapper"],
        mdata.mod["sc"].obs.labels,
    )
    model.project_genes(
        mdata.mod["sc"], mdata.mod["sp"], mdata.mod["sc"].obsm["mapper"]
    )


def test_tangram_errors():
    mdata = _get_mdata()
    Tangram.setup_mudata(
        mdata,
        density_prior_key="rna_count_based_density",
        modalities=modalities,
    )
    with pytest.raises(ValueError):
        Tangram(mdata, constrained=True, target_count=None)

    with pytest.raises(ValueError):
        Tangram.setup_mudata(
            mdata,
            density_prior_key="bad_prior",
            modalities=modalities,
        )
        Tangram(mdata)
