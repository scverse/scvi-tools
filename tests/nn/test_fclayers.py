from __future__ import annotations

import torch
from torch import nn

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
