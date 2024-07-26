from mudata import MuData

from scvi.data import synthetic_iid
from scvi.external import MethylVI


def test_methylvi():
    adata = synthetic_iid()
    adata.layers["mc"] = adata.X
    adata.layers["cov"] = adata.layers["mc"] + 10

    MethylVI.setup_anndata(adata, mc_layer="mc", cov_layer="cov", batch_key="batch")
    vae = MethylVI(
        adata,
    )
    vae.train(3)
    vae.get_elbo(indices=vae.validation_indices)
    vae.get_normalized_methylation()
    vae.get_latent_representation()
    vae.differential_methylation(groupby="labels", group1="label_1")


def test_methylvi_mudata():
    adata1 = synthetic_iid()
    adata1.layers["mc"] = adata1.X
    adata1.layers["cov"] = adata1.layers["mc"] + 10

    adata2 = synthetic_iid()
    adata2.layers["mc"] = adata2.X
    adata2.layers["cov"] = adata2.layers["mc"] + 10

    mdata = MuData({"mod1": adata1, "mod2": adata2})

    MethylVI.setup_mudata(
        mdata,
        mc_layer="mc",
        cov_layer="cov",
        methylation_modalities=["mod1", "mod2"],
        batch_key="batch",
        covariate_modalities={"batch_key": "mod1"},
    )
    vae = MethylVI(
        mdata,
    )
    vae.train(3)
    vae.get_elbo(indices=vae.validation_indices)
    vae.get_normalized_methylation()
    vae.get_latent_representation()
    vae.differential_methylation(groupby="mod1:labels", group1="label_1")
