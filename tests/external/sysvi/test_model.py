import math

import numpy as np
import pandas as pd
from anndata import AnnData
from numpy.testing import assert_raises
from scipy import sparse

from scvi.external import SysVI


def mock_adata():
    # Make random data
    adata = AnnData(
        sparse.csr_matrix(
            np.exp(
                np.concatenate(
                    [
                        np.random.normal(1, 0.5, (200, 5)),
                        np.random.normal(1.1, 0.00237, (200, 5)),
                        np.random.normal(1.3, 0.35, (200, 5)),
                        np.random.normal(2, 0.111, (200, 5)),
                        np.random.normal(2.2, 0.3, (200, 5)),
                        np.random.normal(2.7, 0.01, (200, 5)),
                        np.random.normal(1, 0.001, (200, 5)),
                        np.random.normal(0.00001, 0.4, (200, 5)),
                        np.random.normal(0.2, 0.91, (200, 5)),
                        np.random.normal(0.1, 0.0234, (200, 5)),
                        np.random.normal(0.00005, 0.1, (200, 5)),
                        np.random.normal(0.05, 0.001, (200, 5)),
                        np.random.normal(0.023, 0.3, (200, 5)),
                        np.random.normal(0.6, 0.13, (200, 5)),
                        np.random.normal(0.9, 0.5, (200, 5)),
                        np.random.normal(1, 0.0001, (200, 5)),
                        np.random.normal(1.5, 0.05, (200, 5)),
                        np.random.normal(2, 0.009, (200, 5)),
                        np.random.normal(1, 0.0001, (200, 5)),
                    ],
                    axis=1,
                )
            )
        ),
        var=pd.DataFrame(index=[str(i) for i in range(95)]),
    )
    adata.obs["covariate_cont"] = list(range(200))
    adata.obs["covariate_cat"] = ["a"] * 50 + ["b"] * 50 + ["c"] * 50 + ["d"] * 50
    adata.obs["covariate_cat_emb"] = ["a"] * 50 + ["b"] * 50 + ["c"] * 50 + ["d"] * 50
    adata.obs["system"] = ["a"] * 100 + ["b"] * 50 + ["c"] * 50

    return adata


def test_model():
    adata0 = mock_adata()

    # Run adata setup with all covariates
    SysVI.setup_anndata(
        adata0,
        batch_key="system",
        categorical_covariate_keys=["covariate_cat"],
        categorical_covariate_embed_keys=["covariate_cat_emb"],
        continuous_covariate_keys=["covariate_cont"],
    )

    # Run adata setup transfer
    # TODO ensure this is actually done correctly, not just that it runs through
    adata = mock_adata()
    SysVI.setup_anndata(
        adata,
        batch_key="system",
        categorical_covariate_keys=["covariate_cat"],
        categorical_covariate_embed_keys=["covariate_cat_emb"],
        continuous_covariate_keys=["covariate_cont"],
        covariate_categ_orders=adata0.uns["covariate_categ_orders"],
        covariate_key_orders=adata0.uns["covariate_key_orders"],
        batch_order=adata0.uns["batch_order"],
    )

    # Check that setup of adata without covariates works
    adata_no_cov = mock_adata()
    SysVI.setup_anndata(
        adata_no_cov,
        batch_key="system",
    )
    assert "covariates" not in adata_no_cov.obsm
    assert "covariates_embed" not in adata_no_cov.obsm

    # Model

    # Check that model runs through with standard normal prior
    model = SysVI(adata=adata, prior="standard_normal")
    model.train(max_epochs=2, batch_size=math.ceil(adata.n_obs / 2.0))

    # Check that mode runs through without covariates
    model = SysVI(adata=adata_no_cov)
    model.train(max_epochs=2, batch_size=math.ceil(adata.n_obs / 2.0))

    # Check pre-specifying pseudoinput indices for vamp prior
    _ = SysVI(
        adata=adata,
        prior="vamp",
        pseudoinputs_data_indices=np.array(list(range(5))),
        n_prior_components=5,
    )

    # Check that model runs through with vamp prior without specifying pseudoinput indices,
    # all covariates, and weight scaling
    model = SysVI(adata=adata, prior="vamp")
    model.train(
        max_epochs=2,
        batch_size=math.ceil(adata.n_obs / 2.0),
        log_every_n_steps=1,
        check_val_every_n_epoch=1,
        val_check_interval=1,
        plan_kwargs={
            "log_on_epoch": False,
            "log_on_step": True,
            "loss_weights": {
                "kl_weight": 2,
                "z_distance_cycle_weight": {
                    "weight_start": 1,
                    "weight_end": 3,
                    "point_start": 1,
                    "point_end": 3,
                    "update_on": "step",
                },
            },
        },
    )

    # Embedding

    # Check that embedding default works
    assert (
        model.get_latent_representation(
            adata=adata,
        ).shape[0]
        == adata.shape[0]
    )

    # Ensure that embedding with another adata properly checks if it was setu up correctly
    _ = model.get_latent_representation(adata=adata0)
    with assert_raises(KeyError):
        # TODO could add more check for each property separately
        _ = model.get_latent_representation(adata=adata_no_cov)

    # Check that indices in embedding works
    idx = [1, 2, 3]
    embed = model.get_latent_representation(
        adata=adata,
        indices=idx,
        give_mean=True,
    )
    assert embed.shape[0] == 3

    # Check predicting mean/sample
    np.testing.assert_allclose(
        embed,
        model.get_latent_representation(
            adata=adata,
            indices=idx,
            give_mean=True,
        ),
    )
    with assert_raises(AssertionError):
        np.testing.assert_allclose(
            embed,
            model.get_latent_representation(
                adata=adata,
                indices=idx,
                give_mean=False,
            ),
        )


