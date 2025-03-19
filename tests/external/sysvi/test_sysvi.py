import math

import numpy as np
import pandas as pd
import pytest
from anndata import AnnData
from numpy.testing import assert_raises
from scipy import sparse

from scvi.external import SysVI


def mock_adata():
    """Mock adata for testing."""
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
    adata.obs["batch"] = ["a"] * 100 + ["b"] * 50 + ["c"] * 50

    return adata


@pytest.mark.parametrize(
    (
        "categorical_covariate_keys",
        "continuous_covariate_keys",
        "pseudoinputs_data_indices",
        "embed_categorical_covariates",
        "weight_batches",
    ),
    [
        # Check different covariate combinations
        (["covariate_cat"], ["covariate_cont"], None, False, False),
        (["covariate_cat"], ["covariate_cont"], None, True, False),
        (["covariate_cat"], None, None, False, False),
        (["covariate_cat"], None, None, True, False),
        (None, ["covariate_cont"], None, False, False),
        # Check pre-specifying pseudoinputs
        (None, None, np.array(list(range(5))), False, False),
    ],
)
def test_sysvi_model(
    categorical_covariate_keys,
    continuous_covariate_keys,
    pseudoinputs_data_indices,
    embed_categorical_covariates,
    weight_batches,
):
    """Test model with different input and parameters settings."""
    adata = mock_adata()

    # Run adata setup
    SysVI.setup_anndata(
        adata,
        batch_key="batch",
        categorical_covariate_keys=categorical_covariate_keys,
        continuous_covariate_keys=continuous_covariate_keys,
        weight_batches=weight_batches,
    )

    # Model

    # Check that model runs through with standard normal prior
    model = SysVI(
        adata=adata,
        prior="standard_normal",
        embed_categorical_covariates=embed_categorical_covariates,
    )
    model.train(max_epochs=2, batch_size=math.ceil(adata.n_obs / 2.0))

    # Check that model runs through with vamp prior
    model = SysVI(
        adata=adata,
        prior="vamp",
        pseudoinputs_data_indices=pseudoinputs_data_indices,
        n_prior_components=5,
        embed_categorical_covariates=embed_categorical_covariates,
    )
    model.train(max_epochs=2, batch_size=math.ceil(adata.n_obs / 2.0))

    # Embedding

    # Check that embedding default works
    assert (
        model.get_latent_representation(
            adata=adata,
        ).shape[0]
        == adata.shape[0]
    )


def test_sysvi_latent_representation():
    """Test different parameters for computing later representation."""
    # Train model
    adata = mock_adata()
    SysVI.setup_anndata(
        adata,
        batch_key="batch",
        categorical_covariate_keys=None,
        continuous_covariate_keys=None,
        weight_batches=False,
    )
    model = SysVI(adata=adata, prior="standard_normal")
    model.train(max_epochs=2, batch_size=math.ceil(adata.n_obs / 2.0))

    # Check that specifying indices in embedding works
    idx = [1, 2, 3]
    embed = model.get_latent_representation(
        adata=adata,
        indices=idx,
        give_mean=True,
    )
    assert embed.shape[0] == 3

    # Check predicting mean vs sample
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

    # Test returning distn
    mean, var = model.get_latent_representation(
        adata=adata,
        indices=idx,
        return_dist=True,
    )
    np.testing.assert_allclose(embed, mean)

    model.get_normalized_expression(library_size="latent")
    model.get_normalized_expression(library_size="latent", transform_batch="a")
    model.get_normalized_expression(library_size="latent", indices=[1, 2, 3])


def test_sysvi_warnings():
    """Test that the most important warnings and exceptions are raised."""
    # Train model
    adata = mock_adata()
    SysVI.setup_anndata(
        adata,
        batch_key="batch",
        categorical_covariate_keys=None,
        continuous_covariate_keys=None,
        weight_batches=False,
    )
    model = SysVI(adata=adata, prior="standard_normal")

    # Assert that warning is printed if kl warmup is used
    # Step warmup
    with pytest.warns(Warning) as record:
        model.train(
            max_epochs=2,
            batch_size=math.ceil(adata.n_obs / 2.0),
            plan_kwargs={"n_steps_kl_warmup": 1},
        )
    assert any(
        "The use of KL weight warmup is not recommended in SysVI." in str(rec.message)
        for rec in record
    )
    # Epoch warmup
    with pytest.warns(Warning) as record:
        model.train(
            max_epochs=2,
            batch_size=math.ceil(adata.n_obs / 2.0),
            plan_kwargs={"n_epochs_kl_warmup": 1},
        )
    assert any(
        "The use of KL weight warmup is not recommended in SysVI." in str(rec.message)
        for rec in record
    )

    # Asert that sampling is disabled
    with pytest.raises(NotImplementedError):
        model.module.sample()
