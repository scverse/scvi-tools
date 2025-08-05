import flax.linen as nn
import jax
import jax.numpy as jnp

from scvi.external.mrvi_jax._components import (
    MLP,
    AttentionBlock,
    ConditionalNormalization,
    Dense,
    NormalDistOutputNN,
    ResnetBlock,
)


def test_jaxmrvi_dense():
    key = jax.random.PRNGKey(0)
    x = jnp.ones((20, 10))
    dense = Dense(10)
    params = dense.init(key, x)
    output = dense.apply(params, x)
    assert output.shape == (20, 10)


def test_jaxmrvi_resnetblock():
    key = jax.random.PRNGKey(0)
    x = jnp.ones((20, 10))
    block = ResnetBlock(n_out=30, n_hidden=128, training=True)
    params = block.init(key, x)
    output = block.apply(params, x, mutable=["batch_stats"])
    assert output[0].shape == (20, 30)


def test_jaxmrvi_normalnn():
    key = jax.random.PRNGKey(0)
    key, subkey = jax.random.split(key)
    x = jnp.ones((20, 10))
    nn = NormalDistOutputNN(n_out=30, n_hidden=128, n_layers=3, training=True)
    params = nn.init(key, x)
    output = nn.apply(params, x, mutable=["batch_stats"])
    assert output[0].batch_shape == (20, 30)


def test_jaxmrvi_mlp():
    key = jax.random.PRNGKey(0)
    x = jnp.ones((20, 10))
    mlp = MLP(n_out=30, n_hidden=128, n_layers=3, activation=nn.relu, training=True)
    params = mlp.init(key, x)
    output = mlp.apply(params, x, mutable=["batch_stats"])
    assert output[0].shape == (20, 30)


def test_jaxmrvi_conditionalbatchnorm1d():
    key = jax.random.PRNGKey(0)
    x = jnp.ones((20, 10))
    y = jnp.ones((20, 1))
    conditionalbatchnorm1d = ConditionalNormalization(
        n_features=10, n_conditions=3, normalization_type="batch", training=True
    )
    params = conditionalbatchnorm1d.init(key, x, y)
    output = conditionalbatchnorm1d.apply(params, x, y, mutable=["batch_stats"])
    assert output[0].shape == (20, 10)


def test_jaxmrvi_attention():
    key = jax.random.PRNGKey(0)
    q_vals = jnp.ones((30, 10))
    kv_vals = jnp.ones((30, 10))
    mod = AttentionBlock(query_dim=20, out_dim=40, training=True)
    params = mod.init(key, q_vals, kv_vals)
    mod.apply(params, q_vals, kv_vals, mutable=["batch_stats"])
    q_vals_3d = jnp.ones((3, 30, 10))
    kv_vals_3d = jnp.ones((3, 30, 10))
    r = mod.apply(params, q_vals_3d, kv_vals_3d, mutable=["batch_stats"])
    assert r[0].shape == (3, 30, 40)
