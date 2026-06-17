import math
import os

import numpy as np
import pytest
from anndata import AnnData

import scvi
from scvi.external import DRVI


def mock_adata(n_genes: int = 50, n_batches: int = 2, n_labels: int = 3) -> AnnData:
    """Synthetic count AnnData with batch, labels and extra covariates."""
    adata = scvi.data.synthetic_iid(n_genes=n_genes, n_batches=n_batches, n_labels=n_labels)
    rng = np.random.RandomState(0)
    adata.obs["cont1"] = rng.normal(size=adata.n_obs).astype(np.float32)
    adata.obs["cat1"] = rng.randint(0, 3, size=adata.n_obs).astype(str)
    return adata


@pytest.mark.parametrize(
    ("split_method", "split_aggregation"),
    [
        ("split_map", "logsumexp"),
        ("split_map", "mean"),
        ("split_diag", "logsumexp"),
        ("split_diag", "mean"),
    ],
)
@pytest.mark.parametrize(
    ("categorical_covariate_keys", "continuous_covariate_keys"),
    [
        (None, None),
        (["cat1"], ["cont1"]),
    ],
)
def test_drvi_model(
    split_method,
    split_aggregation,
    categorical_covariate_keys,
    continuous_covariate_keys,
    save_path,
):
    """Train, embed, reconstruct and save under several configurations."""
    adata = mock_adata()
    DRVI.setup_anndata(
        adata,
        batch_key="batch",
        labels_key="labels",
        categorical_covariate_keys=categorical_covariate_keys,
        continuous_covariate_keys=continuous_covariate_keys,
    )

    n_latent = 8
    model = DRVI(
        adata,
        n_latent=n_latent,
        split_method=split_method,
        split_aggregation=split_aggregation,
    )
    model.train(max_epochs=2, batch_size=math.ceil(adata.n_obs / 2.0))

    z = model.get_latent_representation()
    assert z.shape == (adata.n_obs, n_latent)
    # generic VAEMixin / RNASeqMixin methods should work unchanged
    model.get_elbo()
    model.get_reconstruction_error()
    model.get_normalized_expression(n_samples=1)

    dir_path = os.path.join(save_path, "drvi_model/")
    model.save(dir_path, overwrite=True)
    loaded = DRVI.load(dir_path, adata=adata)
    np.testing.assert_allclose(z, loaded.get_latent_representation(), atol=1e-5)


@pytest.mark.parametrize("gene_likelihood", ["nb", "zinb", "poisson"])
def test_drvi_gene_likelihoods(gene_likelihood):
    """The scvi likelihoods are reused unchanged."""
    adata = mock_adata()
    DRVI.setup_anndata(adata, batch_key="batch", labels_key="labels")
    model = DRVI(adata, n_latent=8, gene_likelihood=gene_likelihood)
    model.train(max_epochs=1, batch_size=adata.n_obs)
    assert model.get_latent_representation().shape == (adata.n_obs, 8)


def test_drvi_scarches():
    """scArches transfer learning via the reused ArchesMixin."""
    adata = mock_adata()
    DRVI.setup_anndata(
        adata,
        batch_key="batch",
        labels_key="labels",
        categorical_covariate_keys=["cat1"],
        continuous_covariate_keys=["cont1"],
    )
    model = DRVI(adata, n_latent=8)
    model.train(max_epochs=2, batch_size=math.ceil(adata.n_obs / 2.0))

    query = mock_adata()
    DRVI.prepare_query_anndata(query, model)
    query_model = DRVI.load_query_data(query, model)
    query_model.train(max_epochs=2, batch_size=math.ceil(query.n_obs / 2.0))
    assert query_model.get_latent_representation().shape == (query.n_obs, 8)


@pytest.mark.parametrize(
    ("split_method", "split_aggregation"),
    [("split_map", "logsumexp"), ("split_diag", "mean")],
)
def test_drvi_interpretability(split_method, split_aggregation):
    """Latent-dimension stats and interpretability scores end-to-end."""
    adata = mock_adata()
    DRVI.setup_anndata(adata, batch_key="batch", labels_key="labels")
    n_latent = 8
    model = DRVI(
        adata,
        n_latent=n_latent,
        split_method=split_method,
        split_aggregation=split_aggregation,
    )
    model.train(max_epochs=2, batch_size=adata.n_obs)

    embed = AnnData(model.get_latent_representation(), obs=adata.obs.copy())
    model.set_latent_dimension_stats(embed, adata=adata)
    for col in ["reconstruction_effect", "order", "min", "max", "vanished"]:
        assert col in embed.var

    model.calculate_interpretability_scores(embed, methods="ALL", n_steps=4, n_samples=3)
    assert "OOD_combined_positive" in embed.varm
    assert "IND_max_positive" in embed.varm

    for key in ["OOD_combined", "IND_max"]:
        scores = model.get_interpretability_scores(embed, adata, key=key)
        assert scores.shape[0] == adata.n_vars


def test_drvi_interpretability_requires_full_split():
    """Per-dimension interpretability requires n_split_latent == n_latent."""
    adata = mock_adata()
    DRVI.setup_anndata(adata, batch_key="batch", labels_key="labels")
    model = DRVI(adata, n_latent=8, n_split_latent=4, split_method="split_map")
    model.train(max_epochs=1, batch_size=adata.n_obs)
    embed = AnnData(model.get_latent_representation(), obs=adata.obs.copy())
    with pytest.raises(ValueError, match="one split per latent dimension"):
        model.set_latent_dimension_stats(embed, adata=adata)
