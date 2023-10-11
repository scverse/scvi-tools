import numpy as np
import pytest

from scvi.data import synthetic_iid
from scvi.model import CondSCVI, DestVI


def test_destvi():
    # Step1 learn CondSCVI
    n_latent = 2
    n_labels = 5
    n_layers = 2
    dataset = synthetic_iid(n_labels=n_labels)
    dataset.obs["overclustering_vamp"] = list(range(dataset.n_obs))
    CondSCVI.setup_anndata(dataset, labels_key="labels")
    sc_model = CondSCVI(dataset, n_latent=n_latent, n_layers=n_layers)
    sc_model.train(1, train_size=1)

    # step 2 Check model setup
    DestVI.setup_anndata(dataset, layer=None)

    # Test clustering outside of get_vamp_prior

    # vamp_prior_p>n_largest_cluster to be successful.
    _ = DestVI.from_rna_model(dataset, sc_model, vamp_prior_p=dataset.n_obs)
    # vamp_prior_p<n_largest_cluster leads to value error.
    with pytest.raises(ValueError):
        _ = DestVI.from_rna_model(dataset, sc_model, vamp_prior_p=1)

    del dataset.obs["overclustering_vamp"]

    # step 3 learn destVI with multiple amortization scheme

    for amor_scheme in ["both", "none", "proportion", "latent"]:
        DestVI.setup_anndata(dataset, layer=None)
        # add l1_regularization to cell type proportions
        spatial_model = DestVI.from_rna_model(
            dataset, sc_model, amortization=amor_scheme, l1_reg=50
        )
        spatial_model.view_anndata_setup()
        spatial_model.train(max_epochs=1)
        assert not np.isnan(spatial_model.history["elbo_train"].values[0][0])

        assert spatial_model.get_proportions().shape == (dataset.n_obs, n_labels)
        assert spatial_model.get_gamma(return_numpy=True).shape == (
            dataset.n_obs,
            n_latent,
            n_labels,
        )

        assert spatial_model.get_scale_for_ct("label_0", np.arange(50)).shape == (
            50,
            dataset.n_vars,
        )
