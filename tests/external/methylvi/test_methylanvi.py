from mudata import MuData

from scvi.data import synthetic_iid
from scvi.external import METHYLANVI


def test_methylanvi():
    adata1 = synthetic_iid()
    adata1.layers["mc"] = adata1.X
    adata1.layers["cov"] = adata1.layers["mc"] + 10

    adata2 = synthetic_iid()
    adata2.layers["mc"] = adata2.X
    adata2.layers["cov"] = adata2.layers["mc"] + 10

    mdata = MuData({"mod1": adata1, "mod2": adata2})

    METHYLANVI.setup_mudata(
        mdata,
        mc_layer="mc",
        cov_layer="cov",
        labels_key="labels",
        unlabeled_category="unknown",
        methylation_contexts=["mod1", "mod2"],
        batch_key="batch",
        modalities={"batch_key": "mod1"},
    )
    vae = METHYLANVI(
        mdata,
    )
    vae.train(3)
    vae.get_elbo(indices=vae.validation_indices)
    vae.get_latent_representation()
    vae.predict()
