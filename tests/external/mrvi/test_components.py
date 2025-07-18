import torch

from scvi.external.mrvi._components import (
    AttentionBlock,
    ConditionalNormalization,
    Dense,
    NormalDistOutputNN,
    ResnetBlock,
)


def test_dense():
    torch.manual_seed(0)
    x = torch.ones((20, 10))
    dense = Dense(10, 10)
    output = dense(x)
    assert output.shape == torch.Size([20, 10])


def test_resnetblock():
    torch.manual_seed(0)
    x = torch.ones((20, 10))
    block = ResnetBlock(10, 30, training=True)
    params = block(x)
    assert params.shape == torch.Size([20, 30])


def test_normalnn():
    torch.manual_seed(0)
    x = torch.ones((20, 10))
    nn = NormalDistOutputNN(10, 30, 3, training=True)
    params = nn(x)
    assert params.loc.shape == torch.Size([20, 30])


def test_conditionalbatchnorm1d():
    torch.manual_seed(0)
    x = torch.ones((20, 10))
    y = torch.ones((20, 1))
    conditionalbatchnorm1d = ConditionalNormalization(
        10, 3, normalization_type="batch", training=True
    )
    params = conditionalbatchnorm1d(x, y)
    assert params.shape == torch.Size([20, 10])


def test_attention():
    torch.manual_seed(0)
    q_vals = torch.ones((30, 10))
    kv_vals = torch.ones((30, 10))
    mod = AttentionBlock(query_dim=20, kv_dim=10, out_dim=40, training=True)
    params = mod(q_vals, kv_vals)
    mod.apply(params, q_vals, kv_vals, mutable=["batch_stats"])
    q_vals_3d = torch.ones((3, 30, 10))
    kv_vals_3d = torch.ones((3, 30, 10))
    r = mod(q_vals_3d, kv_vals_3d)
    # r = mod.apply(params, q_vals_3d, kv_vals_3d, mutable=["batch_stats"])
    assert r[0].shape == (3, 30, 40)
