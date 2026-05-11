import pytest
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
        modalities={"batch_key": "mod1", "labels_key": "mod1"},
    )
    vae = METHYLANVI(
        mdata,
    )
    vae.train(3)
    vae.get_elbo(indices=vae.validation_indices)
    vae.get_normalized_methylation()  # Retrieve methylation for all contexts
    vae.get_normalized_methylation(context="mod1")  # Retrieve for specific context
    with pytest.raises(ValueError):  # Should fail when invalid context selected
        vae.get_normalized_methylation(context="mod3")
    vae.get_latent_representation()
    vae.differential_methylation(groupby="mod1:labels", group1="label_1")
    vae.predict()

    predictions, attributions = vae.predict(ig_interpretability=True)
    # let's see an avg of score of top 5 markers for all samples put together
    ig_top_features = attributions.head(5)
    print(ig_top_features)
