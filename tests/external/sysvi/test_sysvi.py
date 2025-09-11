import math
import os

import numpy as np
import pandas as pd
import pytest
from anndata import AnnData
from numpy.testing import assert_raises
from scipy import sparse

from scvi.external import SysVI


def mock_adata(cells_ratio: float = 1):
    """Mock adata for testing."""
    n_cells_base = 200
    n_cells = int(n_cells_base * cells_ratio)
    adata = AnnData(
        sparse.csr_matrix(
            np.exp(
                np.concatenate(
                    [
                        np.random.normal(1, 0.5, (n_cells, 5)),
                        np.random.normal(1.1, 0.00237, (n_cells, 5)),
                        np.random.normal(1.3, 0.35, (n_cells, 5)),
                        np.random.normal(2, 0.111, (n_cells, 5)),
                        np.random.normal(2.2, 0.3, (n_cells, 5)),
                        np.random.normal(2.7, 0.01, (n_cells, 5)),
                        np.random.normal(1, 0.001, (n_cells, 5)),
                        np.random.normal(0.00001, 0.4, (n_cells, 5)),
                        np.random.normal(0.2, 0.91, (n_cells, 5)),
                        np.random.normal(0.1, 0.0234, (n_cells, 5)),
                        np.random.normal(0.00005, 0.1, (n_cells, 5)),
                        np.random.normal(0.05, 0.001, (n_cells, 5)),
                        np.random.normal(0.023, 0.3, (n_cells, 5)),
                        np.random.normal(0.6, 0.13, (n_cells, 5)),
                        np.random.normal(0.9, 0.5, (n_cells, 5)),
                        np.random.normal(1, 0.0001, (n_cells, 5)),
                        np.random.normal(1.5, 0.05, (n_cells, 5)),
                        np.random.normal(2, 0.009, (n_cells, 5)),
                        np.random.normal(1, 0.0001, (n_cells, 5)),
                    ],
                    axis=1,
                )
            )
        ),
        var=pd.DataFrame(index=[str(i) for i in range(95)]),
    )
    adata.obs["covariate_cont"] = list(range(n_cells))
    adata.obs["covariate_cat"] = (
        ["a"] * int(n_cells / 4)
        + ["b"] * int(n_cells / 4)
        + ["c"] * int(n_cells / 4)
        + ["d"] * (n_cells - int(n_cells * 0.75))
    )
    adata.obs["batch"] = (
        ["a"] * int(n_cells / 2)
        + ["b"] * int(n_cells / 4)
        + ["c"] * (n_cells - int(n_cells * 0.75))
    )

    return adata


@pytest.mark.parametrize(
    (
        "prior",
        "categorical_covariate_keys",
        "continuous_covariate_keys",
        "pseudoinputs_data_indices",
        "embed_categorical_covariates",
        "weight_batches",
    ),
    [
        # Check different covariate combinations
        ("vamp", ["covariate_cat"], ["covariate_cont"], None, False, False),
        ("vamp", ["covariate_cat"], ["covariate_cont"], None, True, False),
        ("vamp", ["covariate_cat"], None, None, False, False),
        ("vamp", ["covariate_cat"], None, None, True, False),
        ("vamp", None, ["covariate_cont"], None, False, False),
        # Check alternative priors
        ("standard_normal", ["covariate_cat"], ["covariate_cont"], None, False, False),
        # Check pre-specifying pseudoinputs
        ("vamp", None, None, np.array(list(range(5))), False, False),
        # Check batch weighting
        ("vamp", None, None, None, False, True),
    ],
)
def test_sysvi_model(
    prior,
    categorical_covariate_keys,
    continuous_covariate_keys,
    pseudoinputs_data_indices,
    embed_categorical_covariates,
    weight_batches,
    save_path,
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
    # Check that model runs through
    model = SysVI(
        adata=adata,
        prior=prior,
        pseudoinputs_data_indices=pseudoinputs_data_indices,
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

    # Save model
    dir_path = os.path.join(save_path, "saved_model/")
    model.save(dir_path, overwrite=True)


@pytest.mark.parametrize(
    (
        "prior",
        "embed_categorical_covariates",
    ),
    [
        # Check different covariate representations
        ("vamp", False),
        ("vamp", True),
        # Check different priors
        ("standard_normal", False),
    ],
)
def test_sysvi_scarches(prior, embed_categorical_covariates, save_path):
    # reference adata
    adata = mock_adata()
    SysVI.setup_anndata(
        adata,
        batch_key="batch",
        categorical_covariate_keys=["covariate_cat"],
        continuous_covariate_keys=["covariate_cont"],
    )

    # Reference model
    model = SysVI(
        adata=adata,
        prior=prior,
        embed_categorical_covariates=embed_categorical_covariates,
    )
    model.train(max_epochs=2, batch_size=math.ceil(adata.n_obs / 2.0))

    # Query adata
    adata2 = mock_adata(cells_ratio=0.2)
    # Make it different from reference adata
    adata2.obs["covariate_cat"] = adata2.obs["covariate_cat"].replace({"b": "y", "c": "x"})
    adata2 = adata2[:, np.random.permutation(adata2.var_names)[: adata2.shape[1] - 5]].copy()

    # Make query adata and model
    SysVI.prepare_query_anndata(adata2, model)
    np.testing.assert_equal(adata2.var_names.to_numpy(), adata.var_names.to_numpy())
    model2 = SysVI.load_query_data(adata2, model)

    # Train query model
    model2.train(max_epochs=2, batch_size=math.ceil(adata.n_obs / 2.0))

    # Check that embedding default works
    assert (
        model2.get_latent_representation(
            adata=adata,
        ).shape[0]
        == adata.shape[0]
    )
    assert (
        model2.get_latent_representation(
            adata=adata2,
        ).shape[0]
        == adata2.shape[0]
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
    with pytest.warns(Warning, match="The use of KL weight warmup") as record:
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
    with pytest.warns(Warning, match="The use of KL weight warmup") as record:
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


def test_sysvi_scarches_errors():
    # reference adata
    adata = mock_adata()
    SysVI.setup_anndata(
        adata,
        batch_key="batch",
        categorical_covariate_keys=["covariate_cat"],
        continuous_covariate_keys=["covariate_cont"],
    )

    # Reference model
    model = SysVI(
        adata=adata,
    )
    model.train(max_epochs=2, batch_size=math.ceil(adata.n_obs / 2.0))

    # Query adata
    adata2 = mock_adata()
    # Make it different from reference adata
    adata2.obs["batch"] = adata2.obs["batch"].replace(
        {
            "c": "y",
        }
    )

    # Make query adata and model
    SysVI.prepare_query_anndata(adata2, model)
    with pytest.warns(
        Warning, match="The setting of transfer_batch is disabled in SysVI"
    ) as record:
        with pytest.raises(
            ValueError,
            match="This model does not allow for query having batch categories "
            "missing from the reference.",
        ):
            SysVI.load_query_data(adata2, model, transfer_batch=True)
    assert any(
        "The setting of transfer_batch is disabled in SysVI and is automatically set to False."
        in str(rec.message)
        for rec in record
    )
