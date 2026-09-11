from __future__ import annotations

import pytest
import torch
from torch import nn

import scvi.nn._base_components as base_components_module
from scvi.nn import FCLayers

# ---------------------------------------------------------------------------
# Basic forward behaviour
# ---------------------------------------------------------------------------


def test_forward_2d_output_shape():
    fc = FCLayers(n_in=10, n_out=5, n_layers=2, n_hidden=20)
    x = torch.randn(8, 10)
    assert fc(x).shape == (8, 5)


def test_forward_with_categorical_covariates():
    n_cats = 4
    fc = FCLayers(n_in=10, n_out=5, n_cat_list=[n_cats], use_batch_norm=False, dropout_rate=0.0)
    x = torch.randn(8, 10)
    cat = torch.randint(0, n_cats, (8, 1))
    assert fc(x, cat).shape == (8, 5)


def test_forward_3d_input_with_batch_norm():
    """3D input (n_samples, n_obs, n_in) is handled correctly by the slice-and-cat path."""
    fc = FCLayers(n_in=10, n_out=5, use_batch_norm=True, dropout_rate=0.0)
    fc.eval()
    x = torch.randn(3, 8, 10)
    assert fc(x).shape == (3, 8, 5)


# ---------------------------------------------------------------------------
# _apply_batch_norm: 3D slicing is equivalent to per-sample application
# ---------------------------------------------------------------------------


def test_apply_batch_norm_3d_matches_manual_slicing():
    n_features = 16
    bn = nn.BatchNorm1d(n_features)
    bn.eval()

    fc = FCLayers(n_in=n_features, n_out=n_features, use_batch_norm=True, dropout_rate=0.0)

    x = torch.randn(4, 8, n_features)
    out = fc._apply_batch_norm(bn, x)

    expected = torch.cat([bn(x[i]).unsqueeze(0) for i in range(x.size(0))], dim=0)
    assert out.shape == (4, 8, n_features)
    assert torch.allclose(out, expected)


def test_apply_batch_norm_2d_passthrough():
    n_features = 16
    bn = nn.BatchNorm1d(n_features)
    bn.eval()

    fc = FCLayers(n_in=n_features, n_out=n_features, use_batch_norm=True, dropout_rate=0.0)
    x = torch.randn(8, n_features)
    out = fc._apply_batch_norm(bn, x)

    assert out.shape == (8, n_features)
    assert torch.allclose(out, bn(x))


def test_batchnorm_slice_probe_returns_bool():
    """The probe is cached and answers with a bool on any machine, mps or not."""
    result = base_components_module._mps_supports_batchnorm_slice()
    assert isinstance(result, bool)
    if not torch.backends.mps.is_available():
        assert result is False


@pytest.mark.skipif(not torch.backends.mps.is_available(), reason="requires an MPS device")
def test_apply_batch_norm_3d_on_mps_matches_manual_slicing():
    n_features = 16
    bn = nn.BatchNorm1d(n_features).to("mps")
    bn.eval()

    fc = FCLayers(n_in=n_features, n_out=n_features, use_batch_norm=True, dropout_rate=0.0).to(
        "mps"
    )
    x = torch.randn(4, 8, n_features, device="mps")
    out = fc._apply_batch_norm(bn, x)

    expected = torch.cat([bn(x[i]).unsqueeze(0) for i in range(x.size(0))], dim=0)
    assert out.shape == (4, 8, n_features)
    assert torch.allclose(out, expected)


@pytest.mark.skipif(not torch.backends.mps.is_available(), reason="requires an MPS device")
def test_apply_batch_norm_3d_clones_when_mps_lacks_slicing_support(monkeypatch):
    """Forces the legacy clone-based fallback and checks it still gives the right answer."""
    monkeypatch.setattr(base_components_module, "_mps_supports_batchnorm_slice", lambda: False)

    n_features = 16
    bn = nn.BatchNorm1d(n_features).to("mps")
    bn.eval()

    fc = FCLayers(n_in=n_features, n_out=n_features, use_batch_norm=True, dropout_rate=0.0).to(
        "mps"
    )
    x = torch.randn(4, 8, n_features, device="mps")
    out = fc._apply_batch_norm(bn, x)

    expected = torch.cat([bn(x[i].clone()).unsqueeze(0) for i in range(x.size(0))], dim=0)
    assert out.shape == (4, 8, n_features)
    assert torch.allclose(out, expected)


# ---------------------------------------------------------------------------
# Subclassing seam: _build_layer
# ---------------------------------------------------------------------------


def test_subclass_build_layer_replaces_layer_factory():
    """Overriding _build_layer lets a subclass insert its own layer structure."""

    class CustomFCLayers(FCLayers):
        def _build_layer(self, n_in: int, n_out: int, layer_num: int) -> nn.Sequential:
            return nn.Sequential(
                nn.Linear(n_in + self.n_cov * self.inject_into_layer(layer_num), n_out),
                nn.Tanh(),
            )

    # use_batch_norm and dropout_rate are ignored by the custom factory
    fc = CustomFCLayers(n_in=10, n_out=5, use_batch_norm=True, dropout_rate=0.3)
    layer_types = [type(l) for l in fc.fc_layers[0] if l is not None]

    assert nn.BatchNorm1d not in layer_types
    assert nn.Dropout not in layer_types
    assert nn.Tanh in layer_types

    assert fc(torch.randn(8, 10)).shape == (8, 5)


def test_subclass_build_layer_receives_correct_dims_with_covariates():
    """_build_layer is called with dimensions that already account for covariates on layer 0."""

    recorded = {}

    class CaptureFCLayers(FCLayers):
        def _build_layer(self, n_in: int, n_out: int, layer_num: int) -> nn.Sequential:
            recorded[layer_num] = (n_in, n_out)
            return super()._build_layer(n_in, n_out, layer_num)

    n_cats = 3
    CaptureFCLayers(
        n_in=10,
        n_out=5,
        n_cat_list=[n_cats],
        n_layers=2,
        n_hidden=20,
        use_batch_norm=False,
        dropout_rate=0.0,
    )

    # layer 0: n_in passed is the raw n_in; _build_layer adds n_cov internally
    assert recorded[0] == (10, 20)
    assert recorded[1] == (20, 5)


# ---------------------------------------------------------------------------
# Subclassing seam: _is_linear_layer
# ---------------------------------------------------------------------------


def test_subclass_is_linear_layer_gates_covariate_injection():
    """When _is_linear_layer returns False for a layer, covariates are NOT injected."""

    class NonInjectingFCLayers(FCLayers):
        def _is_linear_layer(self, layer: nn.Module) -> bool:
            return False  # never inject

        def _build_layer(self, n_in: int, n_out: int, layer_num: int) -> nn.Sequential:
            # build without enlarging n_in for covariates so sizes are consistent
            return nn.Sequential(nn.Linear(n_in, n_out))

    n_cats = 4
    fc = NonInjectingFCLayers(
        n_in=10, n_out=5, n_cat_list=[n_cats], use_batch_norm=False, dropout_rate=0.0
    )
    x = torch.randn(8, 10)
    cat = torch.randint(0, n_cats, (8, 1))
    # should not raise: cat tensor is ignored because _is_linear_layer → False
    assert fc(x, cat).shape == (8, 5)


def test_subclass_is_linear_layer_used_by_set_online_update_hooks():
    """set_online_update_hooks registers hooks only on layers identified by _is_linear_layer."""

    class CustomLinear(nn.Linear):
        pass

    class CustomFCLayers(FCLayers):
        def _is_linear_layer(self, layer: nn.Module) -> bool:
            return isinstance(layer, CustomLinear)

        def _build_layer(self, n_in: int, n_out: int, layer_num: int) -> nn.Sequential:
            return nn.Sequential(
                CustomLinear(n_in + self.n_cov * self.inject_into_layer(layer_num), n_out)
            )

    n_cats = 3
    fc = CustomFCLayers(
        n_in=10, n_out=5, n_cat_list=[n_cats], use_batch_norm=False, dropout_rate=0.0
    )
    fc.set_online_update_hooks()
    assert len(fc.hooks) > 0


# ---------------------------------------------------------------------------
# Gradient hook correctness (tests the [... , ] fix for 3-D grad tensors)
# ---------------------------------------------------------------------------


def test_gradient_hook_preserves_categorical_grad_only():
    """After set_online_update_hooks, only categorical weight gradients are non-zero."""
    n_cats = 3
    fc = FCLayers(
        n_in=10,
        n_out=5,
        n_cat_list=[n_cats],
        n_layers=1,
        use_batch_norm=False,
        use_activation=False,
        dropout_rate=0.0,
    )
    fc.set_online_update_hooks()

    x = torch.randn(8, 10)
    cat = torch.randint(0, n_cats, (8, 1))
    fc(x, cat).sum().backward()

    linear = fc.fc_layers[0][0]  # first (only) layer, first sub-module
    grad = linear.weight.grad  # shape: (n_out, n_in + n_cats) = (5, 13)

    # non-categorical columns should be zeroed out by the hook
    assert torch.all(grad[:, :-n_cats] == 0), "non-categorical weight grad should be zero"
    # categorical columns should have non-zero grad (with high probability)
    assert not torch.all(grad[:, -n_cats:] == 0), "categorical weight grad should be non-zero"


# ---------------------------------------------------------------------------
# Subclassing seam: _build_cov_list
# ---------------------------------------------------------------------------


def test_build_cov_list_orders_continuous_before_one_hot_categoricals():
    n_cats = 4
    n_cont = 3
    fc = FCLayers(n_in=10, n_out=5, n_cat_list=[n_cats], n_cont=n_cont)

    cont = torch.randn(8, n_cont)
    cat = torch.randint(0, n_cats, (8, 1))
    cov_list = fc._build_cov_list((cat,), cont)

    assert len(cov_list) == 2
    # continuous covariates come first, categoricals are appended one-hot encoded
    assert torch.allclose(cov_list[0], cont)
    assert cov_list[1].shape == (8, n_cats)
    assert torch.equal(cov_list[1], nn.functional.one_hot(cat.squeeze(-1), n_cats))


def test_build_cov_list_without_covariates_is_empty():
    fc = FCLayers(n_in=10, n_out=5)
    assert fc._build_cov_list((), None) == []


def test_build_cov_list_passes_through_already_one_hot_categoricals():
    """A cat tensor already of width n_cat is used as-is rather than re-encoded."""
    n_cats = 4
    fc = FCLayers(n_in=10, n_out=5, n_cat_list=[n_cats])

    one_hot = nn.functional.one_hot(torch.randint(0, n_cats, (8,)), n_cats)
    cov_list = fc._build_cov_list((one_hot,), None)

    assert len(cov_list) == 1
    assert cov_list[0] is one_hot


def test_build_cov_list_ignores_single_category_covariates():
    """n_cat = 1 carries no information and is dropped from the covariate list."""
    fc = FCLayers(n_in=10, n_out=5, n_cat_list=[1, 3])

    cat_single = torch.zeros(8, 1, dtype=torch.long)
    cat_multi = torch.randint(0, 3, (8, 1))
    cov_list = fc._build_cov_list((cat_single, cat_multi), None)

    assert len(cov_list) == 1
    assert cov_list[0].shape == (8, 3)


def test_build_cov_list_raises_on_too_few_categorical_args():
    fc = FCLayers(n_in=10, n_out=5, n_cat_list=[3, 4])
    with pytest.raises(ValueError, match="doesn't match init"):
        fc._build_cov_list((torch.randint(0, 3, (8, 1)),), None)


def test_build_cov_list_raises_on_missing_categorical_arg():
    fc = FCLayers(n_in=10, n_out=5, n_cat_list=[3])
    with pytest.raises(ValueError, match="cat not provided"):
        fc._build_cov_list((None,), None)


def test_subclass_build_cov_list_overrides_injected_covariates():
    """forward() sources its covariates from _build_cov_list, so a subclass can rewrite them."""

    class ConstantCovFCLayers(FCLayers):
        def _build_cov_list(self, cat_list, cont):
            # ignore the caller's covariates entirely and inject a fixed one-hot instead
            n_obs = cat_list[0].size(0)
            return [nn.functional.one_hot(torch.zeros(n_obs, dtype=torch.long), 3).float()]

    n_cats = 3
    fc = ConstantCovFCLayers(
        n_in=10, n_out=5, n_cat_list=[n_cats], use_batch_norm=False, dropout_rate=0.0
    )
    x = torch.randn(8, 10)

    # two different cat arguments give the same output because they are both discarded
    out_a = fc(x, torch.zeros(8, 1, dtype=torch.long))
    out_b = fc(x, torch.full((8, 1), 2, dtype=torch.long))
    assert torch.allclose(out_a, out_b)


# ---------------------------------------------------------------------------
# Forward kwargs threading to the _apply_layer / _apply_batch_norm seams
# ---------------------------------------------------------------------------


def test_forward_accepts_and_ignores_extra_kwargs():
    """The base layers ignore unknown per-call context instead of raising."""
    fc = FCLayers(n_in=10, n_out=5, n_layers=2, n_hidden=20, dropout_rate=0.0)
    fc.eval()
    x = torch.randn(8, 10)

    assert torch.allclose(fc(x), fc(x, some_context="ignored", other=3))


def test_forward_kwargs_reach_apply_layer_and_apply_batch_norm():
    """Per-call context passed to forward() is threaded to both subclass hooks."""
    seen_layer = []
    seen_batch_norm = []

    class ContextFCLayers(FCLayers):
        def _apply_layer(self, layer, x, cov_list, layer_index, **kwargs):
            seen_layer.append(kwargs)
            return super()._apply_layer(layer, x, cov_list, layer_index, **kwargs)

        def _apply_batch_norm(self, layer, x, **kwargs):
            seen_batch_norm.append(kwargs)
            return super()._apply_batch_norm(layer, x, **kwargs)

    fc = ContextFCLayers(
        n_in=10, n_out=5, n_layers=2, n_hidden=20, use_batch_norm=True, dropout_rate=0.0
    )
    fc.eval()
    fc(torch.randn(8, 10), mode="decode")

    assert len(seen_layer) > 0
    assert all(kw == {"mode": "decode"} for kw in seen_layer)
    assert len(seen_batch_norm) > 0
    assert all(kw == {"mode": "decode"} for kw in seen_batch_norm)


def test_forward_kwargs_can_change_the_output():
    """A subclass can genuinely branch on the per-call context."""

    class ScalingFCLayers(FCLayers):
        def _apply_layer(self, layer, x, cov_list, layer_index, scale=1.0, **kwargs):
            return scale * super()._apply_layer(layer, x, cov_list, layer_index, **kwargs)

    fc = ScalingFCLayers(
        n_in=10,
        n_out=5,
        n_layers=1,
        use_batch_norm=False,
        use_activation=False,
        dropout_rate=0.0,
    )
    fc.eval()
    x = torch.randn(8, 10)

    assert torch.allclose(fc(x, scale=2.0), 2.0 * fc(x, scale=1.0))


def test_forward_without_kwargs_passes_nothing_to_the_hooks():
    """Not passing context leaves the hooks with empty kwargs (no injected defaults)."""
    seen = []

    class ContextFCLayers(FCLayers):
        def _apply_layer(self, layer, x, cov_list, layer_index, **kwargs):
            seen.append(kwargs)
            return super()._apply_layer(layer, x, cov_list, layer_index, **kwargs)

    fc = ContextFCLayers(n_in=10, n_out=5, use_batch_norm=False, dropout_rate=0.0)
    fc(torch.randn(8, 10))

    assert len(seen) > 0
    assert all(kw == {} for kw in seen)


# ---------------------------------------------------------------------------
# Residual skip connections
# ---------------------------------------------------------------------------


def _sequential_blocks(fc: FCLayers) -> list[nn.Sequential]:
    """Rebuild each fc_layers block as a plain Sequential.

    Valid as an independent reference only for 2D input without covariates, where ``_apply_layer``
    and ``_apply_batch_norm`` both reduce to ``layer(x)``.
    """
    return [nn.Sequential(*[l for l in layers if l is not None]) for layers in fc.fc_layers]


def _reference_forward(fc: FCLayers, x: torch.Tensor, residual: bool) -> torch.Tensor:
    h = x
    for i, block in enumerate(_sequential_blocks(fc)):
        out = block(h)
        if residual and i > 0 and out.shape == h.shape:
            out = out + h
        h = out
    return h


def test_residual_defaults_to_off():
    fc = FCLayers(n_in=10, n_out=5)
    assert fc.residual is False


def test_residual_matches_manual_skip_connections():
    """residual=True adds the block input back on every same-width hidden block."""
    torch.manual_seed(0)
    n_hidden = 20
    fc = FCLayers(
        n_in=10,
        n_out=n_hidden,
        n_layers=4,
        n_hidden=n_hidden,
        use_batch_norm=True,
        dropout_rate=0.0,
        residual=True,
    )
    fc.eval()
    x = torch.randn(8, 10)

    assert torch.allclose(fc(x), _reference_forward(fc, x, residual=True), atol=1e-6)


def test_non_residual_matches_manual_forward_without_skips():
    """Sanity check on the reference implementation: residual=False is the plain chain."""
    torch.manual_seed(0)
    n_hidden = 20
    fc = FCLayers(
        n_in=10,
        n_out=n_hidden,
        n_layers=4,
        n_hidden=n_hidden,
        use_batch_norm=True,
        dropout_rate=0.0,
        residual=False,
    )
    fc.eval()
    x = torch.randn(8, 10)

    assert torch.allclose(fc(x), _reference_forward(fc, x, residual=False), atol=1e-6)


def test_residual_changes_the_output():
    """The flag is not a no-op: the same weights give different outputs with and without skips."""
    n_hidden = 20
    kwargs = {
        "n_in": 10,
        "n_out": n_hidden,
        "n_layers": 4,
        "n_hidden": n_hidden,
        "use_batch_norm": True,
        "dropout_rate": 0.0,
    }

    torch.manual_seed(0)
    plain = FCLayers(**kwargs, residual=False).eval()
    torch.manual_seed(0)
    skipped = FCLayers(**kwargs, residual=True).eval()

    x = torch.randn(8, 10)
    assert not torch.allclose(plain(x), skipped(x))


def test_residual_is_a_no_op_for_a_single_layer():
    """Block 0 changes width, so it never gets a skip and n_layers=1 is unaffected."""
    kwargs = {"n_in": 10, "n_out": 5, "n_layers": 1, "dropout_rate": 0.0}

    torch.manual_seed(0)
    plain = FCLayers(**kwargs, residual=False).eval()
    torch.manual_seed(0)
    skipped = FCLayers(**kwargs, residual=True).eval()

    x = torch.randn(8, 10)
    assert torch.allclose(plain(x), skipped(x))


def test_residual_never_skips_the_first_block_even_at_equal_width():
    """The i > 0 guard holds: block 0 gets no skip even when n_in == n_hidden."""
    n_features = 20
    kwargs = {
        "n_in": n_features,
        "n_out": 5,
        "n_layers": 2,
        "n_hidden": n_features,
        "use_batch_norm": True,
        "dropout_rate": 0.0,
    }

    torch.manual_seed(0)
    plain = FCLayers(**kwargs, residual=False).eval()
    torch.manual_seed(0)
    skipped = FCLayers(**kwargs, residual=True).eval()

    x = torch.randn(8, n_features)
    # block 0 maps n_features -> n_features (shapes match) but is excluded by the i > 0 guard,
    # and block 1 maps n_features -> 5, excluded by the shape guard
    assert torch.allclose(plain(x), skipped(x))


def test_residual_skips_blocks_whose_width_changes():
    """With n_hidden != n_out the last block's shapes differ, so no skip is added there."""
    kwargs = {
        "n_in": 10,
        "n_out": 5,
        "n_layers": 2,
        "n_hidden": 20,
        "use_batch_norm": True,
        "dropout_rate": 0.0,
    }

    torch.manual_seed(0)
    plain = FCLayers(**kwargs, residual=False).eval()
    torch.manual_seed(0)
    skipped = FCLayers(**kwargs, residual=True).eval()

    x = torch.randn(8, 10)
    # block 0 is skipped by the i > 0 guard, block 1 by the shape guard -> identical outputs
    assert torch.allclose(plain(x), skipped(x))


def test_residual_with_covariates_and_output_shape():
    n_cats = 4
    n_hidden = 20
    fc = FCLayers(
        n_in=10,
        n_out=n_hidden,
        n_cat_list=[n_cats],
        n_layers=3,
        n_hidden=n_hidden,
        dropout_rate=0.0,
        residual=True,
    )
    fc.eval()
    x = torch.randn(8, 10)
    cat = torch.randint(0, n_cats, (8, 1))

    assert fc(x, cat).shape == (8, n_hidden)


def test_residual_with_3d_input():
    """The skip is added on the (n_samples, n_obs, n_hidden) path too."""
    n_hidden = 20
    fc = FCLayers(
        n_in=10,
        n_out=n_hidden,
        n_layers=3,
        n_hidden=n_hidden,
        use_batch_norm=True,
        dropout_rate=0.0,
        residual=True,
    )
    fc.eval()
    x = torch.randn(3, 8, 10)

    assert fc(x).shape == (3, 8, n_hidden)


def test_residual_gradients_flow_to_the_first_block():
    """The skip path keeps a gradient on early layers even when a later block is zeroed out."""
    n_hidden = 8
    fc = FCLayers(
        n_in=10,
        n_out=n_hidden,
        n_layers=3,
        n_hidden=n_hidden,
        use_batch_norm=False,
        use_activation=False,
        dropout_rate=0.0,
        residual=True,
    )
    # kill the last two blocks: without the skips the output would be identically zero
    for layers in fc.fc_layers[1:]:
        for layer in layers:
            if layer is not None and isinstance(layer, nn.Linear):
                nn.init.zeros_(layer.weight)
                nn.init.zeros_(layer.bias)

    fc(torch.randn(8, 10)).sum().backward()

    first_linear = fc.fc_layers[0][0]
    assert first_linear.weight.grad is not None
    assert not torch.all(first_linear.weight.grad == 0)
