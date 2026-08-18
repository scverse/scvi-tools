"""Tests for VIVS model."""

import numpy as np
import pytest
import torch

from scvi.data import synthetic_iid
from scvi.external.vivs._constants import VIVS_REGISTRY_KEYS


def test_vivs_registry_keys():
    assert VIVS_REGISTRY_KEYS.Y_KEY == "Y"


def test_importance_score_net_mlp_shapes():
    from scvi.external.vivs._module import ImportanceScoreNet

    net = ImportanceScoreNet(n_input=20, n_responses=5, n_hidden=8, loss_type="mse", linear=False)
    x = torch.rand(16, 20)
    y = torch.randn(16, 5)
    out = net(x, y)
    assert out["h"].shape == (16, 5)
    assert out["all_loss"].shape == (16, 5)
    assert out["loss"].dim() == 0


def test_importance_score_net_linear_binary_shapes():
    from scvi.external.vivs._module import ImportanceScoreNet

    net = ImportanceScoreNet(n_input=20, n_responses=1, loss_type="binary", linear=True)
    x = torch.rand(16, 20)
    y = (torch.rand(16, 1) > 0.5).float()
    out = net(x, y)
    assert out["h"].shape == (16, 1)
    assert out["all_loss"].shape == (16, 1)
    assert torch.isfinite(out["loss"])


def test_vivs_module_phase_x_loss():
    from scvi.external.vivs._module import VIVSModule
    from scvi.module.base import LossOutput

    module = VIVSModule(n_input=20, n_responses=5, n_batch=2)
    assert module._phase == "x"
    tensors = {
        "X": torch.rand(8, 20) * 10,
        "batch": torch.randint(0, 2, (8, 1)),
        "labels": torch.zeros(8, 1, dtype=torch.long),
        "Y": torch.randn(8, 5),
    }
    inference_out = module.inference(**module._get_inference_input(tensors))
    generative_out = module.generative(**module._get_generative_input(tensors, inference_out))
    loss_out = module.loss(tensors, inference_out, generative_out)
    assert isinstance(loss_out, LossOutput)
    assert torch.isfinite(loss_out.loss)


def test_vivs_module_phase_xy_loss():
    from scvi.external.vivs._module import VIVSModule

    module = VIVSModule(n_input=20, n_responses=5, n_batch=2)
    module._phase = "xy"
    tensors = {
        "X": torch.rand(8, 20) * 10,
        "batch": torch.randint(0, 2, (8, 1)),
        "labels": torch.zeros(8, 1, dtype=torch.long),
        "Y": torch.randn(8, 5),
    }
    inference_out = module.inference(**module._get_inference_input(tensors))
    generative_out = module.generative(**module._get_generative_input(tensors, inference_out))
    loss_out = module.loss(tensors, inference_out, generative_out)
    assert torch.isfinite(loss_out.loss)
    # xy phase must not touch x_module's parameters
    loss_out.loss.backward()
    for p in module.x_module.parameters():
        assert p.grad is None
    for p in module.xy_module.parameters():
        assert p.grad is not None


@pytest.fixture
def vivs_adata():
    adata = synthetic_iid(n_genes=50)
    return adata


def test_vivs_setup_anndata_and_init(vivs_adata):
    from scvi.external.vivs._model import VIVS

    VIVS.setup_anndata(vivs_adata, y_obsm_key="protein_expression", batch_key="batch")
    model = VIVS(vivs_adata, n_hidden=8, n_latent=4)
    n_proteins = vivs_adata.obsm["protein_expression"].shape[1]
    assert model.module.xy_module.dense_out.out_features == n_proteins
    assert not model.module.x_module_is_pretrained


def test_vivs_train_fresh_vae(vivs_adata):
    from scvi.external.vivs._model import VIVS

    VIVS.setup_anndata(vivs_adata, y_obsm_key="protein_expression", batch_key="batch")
    model = VIVS(vivs_adata, n_hidden=8, n_latent=4)
    model.train(max_epochs=1)
    assert model.is_trained
    assert model.module._phase == "xy"
    assert not any(p.requires_grad for p in model.module.x_module.parameters())


def test_vivs_get_latent_representation(vivs_adata):
    from scvi.external.vivs._model import VIVS

    VIVS.setup_anndata(vivs_adata, y_obsm_key="protein_expression", batch_key="batch")
    model = VIVS(vivs_adata, n_hidden=8, n_latent=4)
    model.train(max_epochs=1)
    z = model.get_latent_representation()
    assert z.shape == (vivs_adata.n_obs, 4)


def test_vivs_pretrained_x_model(vivs_adata):
    from scvi.external.vivs._model import VIVS
    from scvi.model import SCVI

    SCVI.setup_anndata(vivs_adata, batch_key="batch")
    scvi_model = SCVI(vivs_adata, n_hidden=8, n_latent=4)
    scvi_model.train(max_epochs=1)

    VIVS.setup_anndata(vivs_adata, y_obsm_key="protein_expression", batch_key="batch")
    model = VIVS(vivs_adata, x_model=scvi_model)
    assert model.module.x_module_is_pretrained
    assert model.module.x_module is scvi_model.module

    model.train(max_epochs=1)
    assert model.is_trained
    # phase-1 training must have been skipped: x_module's optimizer never ran under VIVS
    assert model.module._phase == "xy"


def test_vivs_untrained_x_model_raises(vivs_adata):
    from scvi.external.vivs._model import VIVS
    from scvi.model import SCVI

    SCVI.setup_anndata(vivs_adata, batch_key="batch")
    scvi_model = SCVI(vivs_adata)  # not trained

    VIVS.setup_anndata(vivs_adata, y_obsm_key="protein_expression", batch_key="batch")
    with pytest.raises(ValueError, match="must already be trained"):
        VIVS(vivs_adata, x_model=scvi_model)


def test_vivs_predict_t(vivs_adata):
    from scvi.external.vivs._model import VIVS

    VIVS.setup_anndata(vivs_adata, y_obsm_key="protein_expression", batch_key="batch")
    model = VIVS(vivs_adata, n_hidden=8, n_latent=4)
    model.train(max_epochs=1)
    t = model.predict_t()
    assert t.shape == (vivs_adata.n_obs, vivs_adata.obsm["protein_expression"].shape[1])
    assert np.isfinite(t).all()


def test_vivs_get_importance_shapes_and_bounds(vivs_adata):
    from scvi.external.vivs._model import VIVS

    VIVS.setup_anndata(vivs_adata, y_obsm_key="protein_expression", batch_key="batch")
    model = VIVS(vivs_adata, n_hidden=8, n_latent=4)
    model.train(max_epochs=1)
    res = model.get_importance(n_mc_samples=10, use_vmap=False)

    n_genes = vivs_adata.n_vars
    n_responses = vivs_adata.obsm["protein_expression"].shape[1]
    assert res["obs_ts"].shape == (n_responses,)
    assert res["null_ts"].shape == (10, n_genes, n_responses)
    assert res["pvalues"].shape == (n_genes, n_responses)
    assert res["padj"].shape == (n_genes, n_responses)
    assert np.all(res["pvalues"] >= 0)
    assert np.all(res["pvalues"] <= 1)
    assert np.all(res["padj"] >= 0)
    assert np.all(res["padj"] <= 1)


def test_vivs_get_importance_vmap_matches_loop(vivs_adata):
    from scvi.external.vivs._model import VIVS

    VIVS.setup_anndata(vivs_adata, y_obsm_key="protein_expression", batch_key="batch")
    model = VIVS(vivs_adata, n_hidden=8, n_latent=4)
    model.train(max_epochs=1)

    torch.manual_seed(0)
    res_loop = model.get_importance(n_mc_samples=5, use_vmap=False)
    torch.manual_seed(0)
    res_vmap = model.get_importance(n_mc_samples=5, use_vmap=True)

    np.testing.assert_allclose(res_loop["obs_ts"], res_vmap["obs_ts"], rtol=1e-4)
    # Null statistics won't match exactly (independent RNG draws per gene under vmap's
    # randomness="different" vs. the loop's shared draw), but should be the same order
    # of magnitude and shape.
    assert res_loop["null_ts"].shape == res_vmap["null_ts"].shape


def test_vivs_get_importance_auto_vmap_threshold(vivs_adata):
    from scvi.external.vivs._model import VIVS

    VIVS.setup_anndata(vivs_adata, y_obsm_key="protein_expression", batch_key="batch")
    model = VIVS(vivs_adata, n_hidden=8, n_latent=4)
    model.train(max_epochs=1)
    # vivs_adata has 50 genes (< 2000), so "auto" must resolve to vmap=True without erroring.
    res = model.get_importance(n_mc_samples=3, use_vmap="auto")
    assert res["pvalues"].shape[0] == vivs_adata.n_vars


def test_vivs_get_cell_scores(vivs_adata):
    from scvi.external.vivs._model import VIVS

    VIVS.setup_anndata(vivs_adata, y_obsm_key="protein_expression", batch_key="batch")
    model = VIVS(vivs_adata, n_hidden=8, n_latent=4)
    model.train(max_epochs=1)
    res = model.get_cell_scores(gene_ids=[0, 1, 2], n_mc_samples=3)
    n_responses = vivs_adata.obsm["protein_expression"].shape[1]
    assert res["tilde_t_mean"].shape[0] == vivs_adata.n_obs
    assert res["obs_t"].shape == (vivs_adata.n_obs, n_responses)
    assert res["obs_t"].shape == res["tilde_t_mean"].shape


def test_vivs_save_load_round_trip(vivs_adata, tmp_path):
    from scvi.external.vivs._model import VIVS

    VIVS.setup_anndata(vivs_adata, y_obsm_key="protein_expression", batch_key="batch")
    model = VIVS(vivs_adata, n_hidden=8, n_latent=4)
    model.train(max_epochs=1)
    preds_before = model.predict_t()

    save_path = str(tmp_path / "vivs_model")
    model.save(save_path, save_anndata=True, overwrite=True)
    loaded_model = VIVS.load(save_path)

    preds_after = loaded_model.predict_t()
    np.testing.assert_allclose(preds_before, preds_after, rtol=1e-4)


def test_select_genes():
    from scvi.external.vivs._utils import select_genes

    adata = synthetic_iid(n_genes=100)
    adata_sub = select_genes(adata, n_top_genes=20)
    assert adata_sub.shape == (adata.n_obs, 20)


def test_select_architecture(vivs_adata):
    from scvi.external.vivs._utils import select_architecture

    best = select_architecture(
        vivs_adata,
        y_obsm_key="protein_expression",
        batch_key="batch",
        xy_model_kwargs_grid={"n_hidden": [4, 8]},
        max_epochs=1,
    )
    assert best["n_hidden"] in (4, 8)


def test_vivs_gene_correlations_and_groupings(vivs_adata):
    from scvi.external.vivs._model import VIVS

    VIVS.setup_anndata(vivs_adata, y_obsm_key="protein_expression", batch_key="batch")
    model = VIVS(vivs_adata, n_hidden=8, n_latent=4)
    model.train(max_epochs=1)

    corr = model.get_gene_correlations()
    n_genes = vivs_adata.n_vars
    assert corr.shape == (n_genes, n_genes)
    assert np.allclose(np.diag(corr), 1.0, atol=1e-2)

    groupings, Z, gene_order = model.get_gene_groupings(n_clusters_list=[5, 10], return_z=True)
    assert len(groupings) == 2
    assert groupings[0].shape == (n_genes,)
    assert len(np.unique(groupings[0])) <= 5
    assert len(gene_order) == n_genes


def test_vivs_get_hier_importance(vivs_adata):
    from scvi.external.vivs._model import VIVS

    VIVS.setup_anndata(vivs_adata, y_obsm_key="protein_expression", batch_key="batch")
    model = VIVS(vivs_adata, n_hidden=8, n_latent=4)
    model.train(max_epochs=1)

    res = model.get_hier_importance(n_clusters_list=[5, 10], batch_size=64)
    assert "pval" in res.data_vars
    assert "padj" in res.data_vars
    assert "cluster_assignment" in res.data_vars
    # 5, 10, plus the always-appended "all genes" (n_genes) resolution
    assert res.sizes["resolution"] == 3
    assert res.sizes["gene_name"] == vivs_adata.n_vars


def test_plot_hier_importance_runs(vivs_adata):
    from scvi.external.vivs._model import VIVS
    from scvi.external.vivs._plotting import plot_hier_importance

    VIVS.setup_anndata(vivs_adata, y_obsm_key="protein_expression", batch_key="batch")
    model = VIVS(vivs_adata, n_hidden=8, n_latent=4)
    model.train(max_epochs=1)
    res = model.get_hier_importance(n_clusters_list=[5, 10], batch_size=64)

    fig = plot_hier_importance(res, plot_fig=True, theme_kwargs={"figure_size": (15, 2)})
    assert fig is not None

    plot_df, labels, breaks = plot_hier_importance(res, plot_fig=False)
    assert isinstance(labels, list)
    assert isinstance(breaks, list)


def test_vivs_import_surface():
    from scvi.external import VIVS
    from scvi.external.vivs import VIVSModule

    assert VIVS.__name__ == "VIVS"
    assert VIVSModule.__name__ == "VIVSModule"
