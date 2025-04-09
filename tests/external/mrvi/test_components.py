import jax
import jax.numpy as jnp

from scvi.external.mrvi._components import (
    AttentionBlock,
    ConditionalNormalization,
    Dense,
    NormalDistOutputNN,
    ResnetBlock,
)


def test_dense():
    key = jax.random.PRNGKey(0)
    x = jnp.ones((20, 10))
    dense = Dense(10)
    params = dense.init(key, x)
    dense.apply(params, x)


def test_resnetblock():
    key = jax.random.PRNGKey(0)
    x = jnp.ones((20, 10))
    block = ResnetBlock(10, 30, training=True)
    params = block.init(key, x)
    block.apply(params, x, mutable=["batch_stats"])


def test_normalnn():
    key = jax.random.PRNGKey(0)
    key, subkey = jax.random.split(key)
    x = jnp.ones((20, 10))
    nn = NormalDistOutputNN(10, 30, 3, training=True)
    params = nn.init(key, x)
    nn.apply(params, x, mutable=["batch_stats"])


def test_conditionalbatchnorm1d():
    key = jax.random.PRNGKey(0)
    x = jnp.ones((20, 10))
    y = jnp.ones((20, 1))
    conditionalbatchnorm1d = ConditionalNormalization(
        10, 3, normalization_type="batch", training=True
    )
    params = conditionalbatchnorm1d.init(key, x, y)
    conditionalbatchnorm1d.apply(params, x, y, mutable=["batch_stats"])


def test_attention():
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
