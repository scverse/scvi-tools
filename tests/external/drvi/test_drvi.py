import math
import os

import numpy as np
import pytest
import torch
from anndata import AnnData
from torch import nn

import scvi
from scvi.external.drvi import DRVI, DRVIModule


def mock_adata(n_genes: int = 50, n_batches: int = 2, n_labels: int = 3) -> AnnData:
    """Synthetic count AnnData with batch, labels and extra covariates."""
    adata = scvi.data.synthetic_iid(n_genes=n_genes, n_batches=n_batches, n_labels=n_labels)
    rng = np.random.RandomState(0)
    adata.obs["cont1"] = rng.normal(size=adata.n_obs).astype(np.float32)
    adata.obs["cat1"] = rng.randint(0, 3, size=adata.n_obs).astype(str)
    return adata


@pytest.mark.parametrize(
    ("split_method", "split_aggregation", "n_split_latent"),
    [
        ("split_map", "logsumexp", None),
        ("split_map", "mean", None),
        ("split_mask", "logsumexp", None),
        ("split_mask", "mean", None),
        # n_split_latent=1: a single split over all latent dims (no disentanglement)
        ("split_map", "logsumexp", 1),
        ("split_mask", "logsumexp", 1),
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
    n_split_latent,
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
        n_split_latent=n_split_latent,
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


@pytest.mark.parametrize(
    "gene_likelihood", ["nb", "pnb", "zinb", "poisson", "normal", "normal_unit_var"]
)
def test_drvi_gene_likelihoods(gene_likelihood):
    """All likelihoods (scvi ones plus DRVI's log-space pnb and normal) train and embed."""
    adata = mock_adata()
    DRVI.setup_anndata(adata, batch_key="batch", labels_key="labels")
    model = DRVI(adata, n_latent=8, gene_likelihood=gene_likelihood)
    model.train(max_epochs=1, batch_size=adata.n_obs)
    assert model.get_latent_representation().shape == (adata.n_obs, 8)
    # reconstruction loss must be finite for every likelihood
    rec = model.get_reconstruction_error()
    assert all(
        torch.isfinite(v).all().item()
        if isinstance(v, torch.Tensor)
        else np.isfinite(np.asarray(v)).all()
        for v in rec.values()
    )


def test_drvi_registry_init():
    """Model can be initialized from a registry (no in-memory AnnData), matching scvi's
    datamodule/annbatch interface."""
    adata = mock_adata()
    DRVI.setup_anndata(
        adata,
        batch_key="batch",
        labels_key="labels",
        categorical_covariate_keys=["cat1"],
        continuous_covariate_keys=["cont1"],
    )
    model = DRVI(adata, n_latent=8)
    registry = model.adata_manager.registry

    model_from_registry = DRVI(registry=registry, n_latent=8)
    assert model_from_registry.module.n_input == model.module.n_input
    assert model_from_registry.module.n_batch == model.module.n_batch
    assert model_from_registry.module.n_latent == 8


@pytest.mark.parametrize("activation_fn", ["elu", "relu"])
def test_drvi_activation_fn(activation_fn):
    """Hidden-layer activation is configurable (default ELU) in both encoder and decoder."""
    expected = {"elu": "ELU", "relu": "ReLU"}[activation_fn]
    adata = mock_adata()
    DRVI.setup_anndata(adata, batch_key="batch", labels_key="labels")
    model = DRVI(adata, n_latent=8, activation_fn=activation_fn)

    enc = {type(m).__name__ for m in model.module.z_encoder.encoder.modules()}
    dec = {type(m).__name__ for m in model.module.decoder.px_decoder.modules()}
    assert expected in enc
    assert expected in dec

    model.train(max_epochs=1, batch_size=adata.n_obs)
    assert model.get_latent_representation().shape == (adata.n_obs, 8)


@pytest.mark.parametrize(
    ("spec", "expected_type", "attr", "value"),
    [
        ("ReLU", nn.ReLU, None, None),
        ("ELU", nn.ELU, "alpha", 1.0),
        ("ELU_0.5", nn.ELU, "alpha", 0.5),
        ("nn.ELU_0.5", nn.ELU, "alpha", 0.5),
    ],
)
def test_drvi_resolve_mean_activation(spec, expected_type, attr, value):
    """The resolver maps a torch.nn class name (with optional ``_<arg>``) to an instance."""
    from scvi.external.drvi._module import _resolve_mean_activation

    act = _resolve_mean_activation(spec)
    assert isinstance(act, expected_type)
    if attr is not None:
        assert getattr(act, attr) == value
    # None and an nn.Module instance round-trip
    assert isinstance(_resolve_mean_activation(None), nn.Identity)
    passthrough = nn.ReLU()
    assert _resolve_mean_activation(passthrough) is passthrough


def test_drvi_mean_activation(save_path):
    """mean_activation wraps the latent mean head only, constrains q_m, and round-trips."""
    adata = mock_adata()
    DRVI.setup_anndata(adata, batch_key="batch", labels_key="labels")
    model = DRVI(adata, n_latent=8, mean_activation="ReLU", use_observed_lib_size=False)

    # wraps the latent mean head (Sequential(Linear, ReLU))
    assert isinstance(model.module.z_encoder.mean_encoder[-1], nn.ReLU)

    model.train(max_epochs=2, batch_size=adata.n_obs)
    # relu on q_m => non-negative latent means (latent_distribution="normal" => identity transform)
    z_mean = model.get_latent_representation(give_mean=True)
    assert (z_mean >= 0).all()

    # the wrapped mean head round-trips through save/load (param names stay consistent)
    dir_path = os.path.join(save_path, "drvi_mean_activation/")
    model.save(dir_path, overwrite=True)
    loaded = DRVI.load(dir_path, adata=adata)
    np.testing.assert_allclose(z_mean, loaded.get_latent_representation(give_mean=True), atol=1e-5)


def test_drvi_batch_embedding():
    """Batch embedding (EmbeddingMixin) is supported: the embedded batch is injected into each
    split, and get_batch_representation works."""
    adata = mock_adata()
    DRVI.setup_anndata(adata, batch_key="batch", labels_key="labels")
    model = DRVI(adata, n_latent=8, batch_representation="embedding")
    model.train(max_epochs=2, batch_size=adata.n_obs)
    assert model.get_latent_representation().shape == (adata.n_obs, 8)
    assert model.get_batch_representation().shape[0] == adata.n_obs


@pytest.mark.parametrize(
    ("decoder_reuse_weights", "n_body_per_split", "head_is_shared"),
    [
        ("everywhere", 0, True),
        ("hidden", 0, False),
        ("last", 2, True),
        ("hidden_except_first", 1, True),
        ("nowhere", 2, False),
    ],
)
def test_drvi_decoder_reuse_weights(decoder_reuse_weights, n_body_per_split, head_is_shared):
    """decoder_reuse_weights selects which decoder layers are shared (nn.Linear) vs per-split
    (StackedLinearLayer), across the FC body (2 layers here) and the parameter heads."""
    from stacked_linear import StackedLinearLayer

    adata = mock_adata()
    DRVI.setup_anndata(adata, batch_key="batch", labels_key="labels")
    model = DRVI(adata, n_latent=8, n_layers=2, decoder_reuse_weights=decoder_reuse_weights)

    decoder = model.module.decoder
    body_per_split = sum(isinstance(m, StackedLinearLayer) for m in decoder.px_decoder.modules())
    assert body_per_split == n_body_per_split
    assert isinstance(decoder.px_scale_decoder, nn.Linear) == head_is_shared

    model.train(max_epochs=1, batch_size=adata.n_obs)
    assert model.get_latent_representation().shape == (adata.n_obs, 8)


def test_drvi_module_decoder_regularization_options_follow_user_settings():
    """Decoder layer norm and dropout follow the user's DRVI module settings."""
    module = DRVIModule(
        n_input=10,
        n_batch=1,
        n_labels=1,
        n_latent=4,
        n_hidden=8,
        n_split_latent=2,
        use_layer_norm="none",
        dropout_rate=0.25,
    )

    decoder_modules = list(module.decoder.px_decoder.modules())

    dropouts = [m for m in decoder_modules if isinstance(m, nn.Dropout)]
    assert dropouts
    assert all(m.p == 0.25 for m in dropouts)


def test_drvi_size_factor_key():
    """size_factor_key is registered (full SCVI setup_anndata) and toggles use_size_factor_key."""
    adata = mock_adata()
    adata.obs["size_factor"] = np.asarray(adata.X.sum(1)).ravel().astype(np.float32)
    DRVI.setup_anndata(
        adata, batch_key="batch", labels_key="labels", size_factor_key="size_factor"
    )
    model = DRVI(adata, n_latent=8)
    assert model.module.use_size_factor_key
    model.train(max_epochs=2, batch_size=adata.n_obs)
    assert model.get_latent_representation().shape == (adata.n_obs, 8)


def test_drvi_module_requires_size_factor_when_configured():
    """A module configured with size_factor_key fails clearly if the tensor is missing."""
    module = DRVIModule(
        n_input=10,
        n_batch=1,
        n_labels=1,
        n_latent=4,
        n_hidden=8,
        n_split_latent=2,
        use_size_factor_key=True,
    )

    with pytest.raises(ValueError, match="use_size_factor_key=True"):
        module.generative(
            z=torch.randn(3, 4),
            library=torch.zeros(3, 1),
            batch_index=torch.zeros(3, 1, dtype=torch.long),
            y=torch.zeros(3, 1, dtype=torch.long),
        )


def test_drvi_minified():
    """Minified mode (BaseMinifiedModeModelClass) works end-to-end."""
    adata = mock_adata()
    DRVI.setup_anndata(adata, batch_key="batch", labels_key="labels")
    model = DRVI(adata, n_latent=8)
    model.train(max_epochs=2, batch_size=adata.n_obs)

    qzm, qzv = model.get_latent_representation(give_mean=True, return_dist=True)
    adata.obsm["X_latent_qzm"] = qzm
    adata.obsm["X_latent_qzv"] = qzv
    model.minify_adata()
    assert model.minified_data_type is not None
    assert model.get_latent_representation().shape == (adata.n_obs, 8)
    assert model.get_normalized_expression().shape == (adata.n_obs, adata.n_vars)


@pytest.mark.parametrize("dispersion", ["gene", "gene-batch", "gene-cell"])
def test_drvi_pnb_dispersion(dispersion):
    """The parametrized NB works across dispersion modes."""
    adata = mock_adata()
    DRVI.setup_anndata(adata, batch_key="batch", labels_key="labels")
    model = DRVI(adata, n_latent=8, gene_likelihood="pnb", dispersion=dispersion)
    model.train(max_epochs=1, batch_size=adata.n_obs)
    assert model.get_latent_representation().shape == (adata.n_obs, 8)


def test_lognegativebinomial_matches_scvi_nb():
    """LogNegativeBinomial.log_prob equals scvi NegativeBinomial.log_prob for matching params."""
    import torch

    from scvi.distributions import NegativeBinomial
    from scvi.external.drvi import LogNegativeBinomial

    torch.manual_seed(0)
    log_m = torch.randn(8, 16) * 3.0
    log_r = torch.randn(8, 16) * 2.0
    value = torch.randint(0, 50, (8, 16)).float()

    lp_log = LogNegativeBinomial(log_m=log_m, log_r=log_r).log_prob(value)
    lp_nb = NegativeBinomial(mu=torch.exp(log_m), theta=torch.exp(log_r)).log_prob(value)
    assert torch.allclose(lp_log, lp_nb, atol=1e-3, rtol=1e-4)

    # stable (finite) even where the linear-space mean under/overflows
    log_m_extreme = torch.tensor([[-60.0, 40.0, 0.0]])
    log_r_extreme = torch.zeros(1, 3)
    val = torch.tensor([[0.0, 5.0, 3.0]])
    assert torch.isfinite(
        LogNegativeBinomial(log_m=log_m_extreme, log_r=log_r_extreme).log_prob(val)
    ).all()

    # reassigning mu (as RNASeqMixin importance weighting does) keeps log_prob consistent
    dist = LogNegativeBinomial(log_m=torch.zeros(1, 3), log_r=torch.zeros(1, 3))
    dist.mu = torch.full((1, 3), 5.0)
    val2 = torch.tensor([[1.0, 2.0, 3.0]])
    expected = NegativeBinomial(mu=dist.mu, theta=torch.ones(1, 3)).log_prob(val2)
    assert torch.allclose(dist.log_prob(val2), expected, atol=1e-4)


def _run_query_to_reference(**model_kwargs):
    """scArches query->reference run; returns (reference_change, query_change).

    Holds out ``batch_0`` as a new query batch, trains the reference model, transfers via
    ``load_query_data`` and finetunes on the query. Then measures how much the reference
    and query latent representations change.
    """
    adata = mock_adata(n_batches=3)
    ref = adata[adata.obs["batch"] != "batch_0"].copy()
    query = adata[adata.obs["batch"] == "batch_0"].copy()

    DRVI.setup_anndata(ref, batch_key="batch", labels_key="labels")
    model = DRVI(ref, n_latent=8, encode_covariates=True, **model_kwargs)
    model.train(max_epochs=2, batch_size=ref.n_obs)
    latent_reference = model.get_latent_representation(ref)

    DRVI.prepare_query_anndata(query, model)
    transfer = DRVI.load_query_data(query, model)
    train_kwargs = {"plan_kwargs": {"lr": 0.1, "weight_decay": 0.0}}
    transfer.train(max_epochs=2, batch_size=query.n_obs, **train_kwargs)
    latent_query = transfer.get_latent_representation(query)
    transfer.train(max_epochs=2, batch_size=query.n_obs, **train_kwargs)

    reference_change = np.sum((latent_reference - transfer.get_latent_representation(ref)) ** 2)
    query_change = np.sum((latent_query - transfer.get_latent_representation(query)) ** 2)
    return reference_change, query_change


@pytest.mark.parametrize(
    "model_kwargs",
    [
        # batch modeling: one-hot (also the default used by the configs below)
        {"batch_representation": "one-hot"},
        # weight-sharing strategies
        {"decoder_reuse_weights": "everywhere"},
        {"decoder_reuse_weights": "hidden"},
        {"decoder_reuse_weights": "last"},
        {"decoder_reuse_weights": "hidden_except_first"},
        {"decoder_reuse_weights": "nowhere"},
        # splitting strategies (split_map is the default of the configs above)
        {"split_method": "split_mask"},
        {"n_split_latent": 4},
    ],
    ids=[
        "one-hot",
        "reuse-everywhere",
        "reuse-hidden",
        "reuse-last",
        "reuse-hidden_except_first",
        "reuse-nowhere",
        "split_mask",
        "n_split_4",
    ],
)
def test_drvi_query_to_reference_mapping(model_kwargs):
    """scArches query->reference mapping keeps reference embeddings intact while the query
    embeddings are updated."""
    reference_change, query_change = _run_query_to_reference(**model_kwargs)
    assert reference_change < 1e-6
    assert query_change > 1e-3


def test_drvi_query_to_reference_mapping_embedding_batch():
    """With ``batch_representation="embedding"``, scArches keeps the reference latent intact
    while updating the query."""
    reference_change, query_change = _run_query_to_reference(batch_representation="embedding")
    assert reference_change < 1e-6  # reference embeddings intact
    assert query_change > 1e-3  # query updated by transfer training


def test_drvi_query_to_reference_freezes_reference_embedding_rows():
    """With ``batch_representation="embedding"``, scArches trains only the new query row."""
    adata = mock_adata(n_batches=3)
    ref = adata[adata.obs["batch"] != "batch_0"].copy()
    query = adata[adata.obs["batch"] == "batch_0"].copy()

    DRVI.setup_anndata(ref, batch_key="batch", labels_key="labels")
    model = DRVI(ref, n_latent=8, encode_covariates=True, batch_representation="embedding")
    model.train(max_epochs=2, batch_size=ref.n_obs)
    n_old = model.summary_stats.n_batch

    DRVI.prepare_query_anndata(query, model)
    transfer = DRVI.load_query_data(query, model)
    emb = transfer.module.get_embedding(scvi.REGISTRY_KEYS.BATCH_KEY)
    assert emb.num_embeddings == n_old + 1
    assert emb.weight.requires_grad

    before = emb.weight.detach().clone()
    transfer.train(
        max_epochs=3, batch_size=query.n_obs, plan_kwargs={"lr": 0.1, "weight_decay": 0.0}
    )
    after = emb.weight.detach()

    assert ((after[:n_old] - before[:n_old]) ** 2).sum().item() < 1e-8
    assert ((after[n_old:] - before[n_old:]) ** 2).sum().item() > 1e-3


@pytest.mark.parametrize("deeply_inject_covariates", [True, False])
def test_drvi_query_to_reference_freezes_per_split_decoder_weights(deeply_inject_covariates):
    """With weight-sharing off (per-split ``StackedLinearLayer``, 3d weights), scArches keeps the
    reference decoder weights frozen while updating only the new covariate (query batch) column.
    """
    from stacked_linear import StackedLinearLayer

    adata = mock_adata(n_batches=3)
    ref = adata[adata.obs["batch"] != "batch_0"].copy()
    query = adata[adata.obs["batch"] == "batch_0"].copy()
    DRVI.setup_anndata(ref, batch_key="batch", labels_key="labels")
    model = DRVI(
        ref,
        n_latent=8,
        n_layers=3,
        decoder_reuse_weights="nowhere",
        deeply_inject_covariates=deeply_inject_covariates,
    )
    model.train(max_epochs=2, batch_size=ref.n_obs)

    DRVI.prepare_query_anndata(query, model)
    transfer = DRVI.load_query_data(query, model)

    # Checking all per-split (StackedLinearLayer) layers
    px = transfer.module.decoder.px_decoder
    stacked = [
        (i, group[0])
        for i, group in enumerate(px.fc_layers)
        if isinstance(group[0], StackedLinearLayer)
    ]
    assert stacked  # at least one per-split layer exists

    # the new query batch is appended last -> its covariate column is the last input column
    before = {i: layer.weight.detach().clone() for i, layer in stacked}
    transfer.train(
        max_epochs=3, batch_size=query.n_obs, plan_kwargs={"lr": 0.1, "weight_decay": 0.0}
    )

    n_updated = 0
    for i, layer in stacked:
        delta = (layer.weight.detach() - before[i]) ** 2
        if px.inject_into_layer(i):
            assert delta[..., :-1].sum().item() < 1e-8  # reference weights frozen
            assert delta[..., -1:].sum().item() > 1e-3  # new query-batch column updated
            n_updated += 1
        else:
            assert delta.sum().item() < 1e-8  # no covariate column -> fully frozen


@pytest.mark.parametrize("batch_representation", ["one-hot", "embedding"])
@pytest.mark.parametrize(
    ("split_method", "split_aggregation"),
    [("split_map", "logsumexp"), ("split_mask", "mean")],
)
def test_drvi_interpretability(split_method, split_aggregation, batch_representation):
    """Latent-dimension stats and interpretability scores end-to-end."""
    adata = mock_adata()
    DRVI.setup_anndata(adata, batch_key="batch", labels_key="labels")
    n_latent = 8
    model = DRVI(
        adata,
        n_latent=n_latent,
        split_method=split_method,
        split_aggregation=split_aggregation,
        batch_representation=batch_representation,
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


@pytest.mark.parametrize("split_aggregation", ["logsumexp", "mean"])
def test_drvi_interpretability_split_smaller_than_latent(split_aggregation):
    """Non-directional IND interpretability works when splits group several latent dims.

    With ``n_split_latent < n_latent`` each split covers a contiguous chunk of latent dims, so the
    per-split effects are ``(n_split, n_genes)`` (i.e. how much genes are affected by each split).
    """
    adata = mock_adata()
    DRVI.setup_anndata(adata, batch_key="batch", labels_key="labels")
    n_latent, n_split = 8, 4
    model = DRVI(
        adata,
        n_latent=n_latent,
        n_split_latent=n_split,
        split_method="split_map",
        split_aggregation=split_aggregation,
    )
    model.train(max_epochs=2, batch_size=adata.n_obs)

    # per-split reconstruction effect: one value per split
    recon = model.get_reconstruction_effect_of_each_split(directional=False)
    assert recon.shape == (n_split,)

    # per-split, per-gene effects: every aggregation is (n_split, n_genes)
    effects = model.get_effect_of_splits_within_distribution(directional=False)
    for agg, value in effects.items():
        assert value.shape == (n_split, adata.n_vars), agg

    # directional in-distribution has no clean per-split direction here -> explicit error
    with pytest.raises(NotImplementedError):
        model.get_effect_of_splits_within_distribution(directional=True)

    # set_latent_dimension_stats broadcasts each split's reconstruction effect to its latent dims
    embed = AnnData(model.get_latent_representation(), obs=adata.obs.copy())
    model.set_latent_dimension_stats(embed, adata=adata)
    recon = embed.var["reconstruction_effect"].to_numpy()
    assert recon.shape == (n_latent,)
    # each split covers n_latent // n_split contiguous dims that share the split's effect
    d = n_latent // n_split
    for i in range(n_split):
        block = recon[i * d : (i + 1) * d]
        assert np.allclose(block, block[0])


@pytest.mark.parametrize("gene_likelihood", ["nb", "pnb", "zinb"])
def test_drvi_rnaseq_mixin(gene_likelihood):
    """The inherited :class:`~scvi.model.base.RNASeqMixin` methods work with DRVI.

    Covers the parametrized log-space NB (``pnb``) too, exercising
    :class:`~scvi.external.drvi.LogNegativeBinomial`'s ``mu`` / ``theta`` / ``scale`` interface and
    the ``mu``-reassignment used by importance weighting.
    """
    n_genes, n_latent = 50, 8
    adata = mock_adata(n_genes=n_genes)
    DRVI.setup_anndata(adata, batch_key="batch", labels_key="labels")
    model = DRVI(adata, n_latent=n_latent, gene_likelihood=gene_likelihood)
    model.train(max_epochs=2, batch_size=adata.n_obs)

    # get_normalized_expression: defaults, multi-sample, per-batch, gene subset, latent library
    norm = model.get_normalized_expression()
    assert norm.shape == (adata.n_obs, n_genes)
    samples = model.get_normalized_expression(n_samples=3, return_mean=False, return_numpy=True)
    assert samples.shape == (3, adata.n_obs, n_genes)
    tb = model.get_normalized_expression(transform_batch=["batch_0"])
    assert tb.shape == (adata.n_obs, n_genes)
    subset = model.get_normalized_expression(gene_list=adata.var_names[:5].tolist())
    assert subset.shape == (adata.n_obs, 5)
    latent_lib = model.get_normalized_expression(library_size="latent")
    assert latent_lib.shape == (adata.n_obs, n_genes)

    # importance-weighted expression (reassigns px.mu before log_prob -> must stay consistent)
    importance = model.get_normalized_expression(
        n_samples=5, weights="importance", return_numpy=True
    )
    assert importance.shape == (adata.n_obs, n_genes)
    assert np.isfinite(importance).all()

    # likelihood parameters
    params = model.get_likelihood_parameters()
    assert "mean" in params
    assert params["mean"].shape == (adata.n_obs, n_genes)

    # posterior predictive sampling
    pps = model.posterior_predictive_sample()
    assert pps.shape == (adata.n_obs, n_genes)

    # observed library size (DRVI has no library encoder, matching SCVI's observed-lib default)
    lib = model.get_latent_library_size(give_mean=False)
    assert lib.shape == (adata.n_obs, 1)

    # feature-feature correlation matrix
    corr = model.get_feature_correlation_matrix(n_samples=3)
    assert corr.shape == (n_genes, n_genes)

    # differential expression (change mode, one-vs-rest over labels)
    de = model.differential_expression(groupby="labels", mode="change", silent=True)
    assert de.shape[0] > 0
    assert "proba_de" in de.columns
