import pytest
import torch
from torch import nn

from scvi.external.mrvi_torch._components import (
    MLP,
    AttentionBlock,
    ConditionalNormalization,
    Dense,
    NormalDistOutputNN,
    ResnetBlock,
)


def _init_weights(m):
    if isinstance(m, nn.Linear):
        nn.init.xavier_uniform_(m.weight)
        if m.bias is not None:
            nn.init.zeros_(m.bias)


def test_torchmrvi_dense():
    torch.manual_seed(0)
    x = torch.ones((20, 10))
    dense = Dense(10, 10)
    output = dense(x)
    assert output.shape == torch.Size([20, 10])


def test_torchmrvi_resnetblock():
    torch.manual_seed(0)
    x = torch.ones((20, 10))
    block = ResnetBlock(n_in=10, n_out=30, n_hidden=128)
    block.apply(_init_weights)
    block.train()
    params = block(x)
    assert params.shape == torch.Size([20, 30])


def test_torchmrvi_normalnn():
    torch.manual_seed(0)
    x = torch.ones((20, 10))
    nn = NormalDistOutputNN(n_in=10, n_out=30, n_hidden=128, n_layers=3)
    params = nn(x)
    assert params.loc.shape == torch.Size([20, 30])


def test_torchmrvi_mlp():
    torch.manual_seed(0)
    x = torch.ones((20, 10))
    mlp = MLP(n_in=10, n_out=30, n_hidden=128, n_layers=3, activation=nn.ReLU())
    mlp.train()
    output = mlp(x)
    assert output.shape == torch.Size([20, 30])


@pytest.mark.parametrize("training", [True, False])
def test_torchmrvi_conditionalbatchnorm1d(training):
    torch.manual_seed(0)
    x = torch.ones((20, 10))
    y = torch.ones((20, 1))
    conditionalbatchnorm1d = ConditionalNormalization(
        n_features=10,
        n_conditions=3,
        normalization_type="batch",
    )
    conditionalbatchnorm1d.train()
    params = conditionalbatchnorm1d(x, y, training=training)
    assert params.shape == torch.Size([20, 10])


def test_torchmrvi_attention():
    torch.manual_seed(0)
    q_vals = torch.ones((30, 10))
    kv_vals = torch.ones((30, 10))
    mod = AttentionBlock(query_dim=10, kv_dim=10, out_dim=40)
    out = mod(q_vals, kv_vals)
    assert out.shape == (30, 40)
    q_vals_3d = torch.ones((3, 30, 10))
    kv_vals_3d = torch.ones((3, 30, 10))
    r = mod(q_vals_3d, kv_vals_3d)
    assert r.shape == (3, 30, 40)
