"""Numerical equivalence test between JAX and PyTorch MRVI implementations.

This test initializes both implementations with identical weights,
feeds them the same inputs, and verifies that outputs match numerically.
"""

from __future__ import annotations

import os

# Force full float32 precision (disable TF32 on Ampere GPUs)
os.environ["JAX_DEFAULT_MATMUL_PRECISION"] = "float32"

import numpy as np
import torch

torch.backends.cuda.matmul.allow_tf32 = False
torch.backends.cudnn.allow_tf32 = False
torch.set_float32_matmul_precision("highest")

import jax  # noqa: E402
import jax.numpy as jnp  # noqa: E402

jax.config.update("jax_default_matmul_precision", "float32")


# For simple components (Dense, MLP, ResnetBlock, etc.)
ATOL = 1e-5
RTOL = 1e-4

# For compound components (AttentionBlock) where softmax numerical differences
# amplify through downstream MLPs in float32
ATOL_COMPOUND = 5e-4
RTOL_COMPOUND = 5e-3


def _transfer_dense_weights(jax_params, torch_module, prefix=""):
    """Transfer Dense/Linear weights from JAX to PyTorch."""
    kernel = np.array(jax_params["kernel"])  # (in, out)
    torch_module.weight.data = torch.tensor(kernel.T)  # PyTorch: (out, in)
    if "bias" in jax_params:
        torch_module.bias.data = torch.tensor(np.array(jax_params["bias"]))


def _transfer_layernorm_weights(jax_params, torch_module):
    """Transfer LayerNorm weights from JAX to PyTorch."""
    if "scale" in jax_params:
        torch_module.weight.data = torch.tensor(np.array(jax_params["scale"]))
    if "bias" in jax_params:
        torch_module.bias.data = torch.tensor(np.array(jax_params["bias"]))


def _transfer_embedding_weights(jax_params, torch_module):
    """Transfer Embedding weights from JAX to PyTorch."""
    torch_module.weight.data = torch.tensor(np.array(jax_params["embedding"]))


def _transfer_resnet_block_weights(jax_params, torch_rb):
    """Transfer ResnetBlock weights from JAX to PyTorch.

    Handles the variable number of Dense layers depending on n_in vs n_hidden.
    """
    _transfer_dense_weights(jax_params["Dense_0"], torch_rb.fc1)
    _transfer_layernorm_weights(jax_params["LayerNorm_0"], torch_rb.layer_norm1)

    if torch_rb.fc_match is not None:
        # n_in != n_hidden: Dense_1 = fc_match, Dense_2 = fc2
        _transfer_dense_weights(jax_params["Dense_1"], torch_rb.fc_match)
        _transfer_dense_weights(jax_params["Dense_2"], torch_rb.fc2)
    else:
        # n_in == n_hidden: Dense_1 = fc2 (skip is identity)
        _transfer_dense_weights(jax_params["Dense_1"], torch_rb.fc2)
    _transfer_layernorm_weights(jax_params["LayerNorm_1"], torch_rb.layer_norm2)


def _transfer_mlp_weights(jax_params, torch_mlp):
    """Transfer MLP weights from JAX to PyTorch."""
    for i in range(len(torch_mlp.resnet_blocks)):
        _transfer_resnet_block_weights(jax_params[f"ResnetBlock_{i}"], torch_mlp.resnet_blocks[i])
    _transfer_dense_weights(jax_params["Dense_0"], torch_mlp.fc)


def test_resnet_block_numerical():
    """Test ResnetBlock numerical equivalence with transferred weights."""
    from scvi.external.mrvi_jax._components import ResnetBlock as JaxResnetBlock
    from scvi.external.mrvi_torch._components import ResnetBlock as TorchResnetBlock

    n_in = 32
    n_out = 16
    n_hidden = 32
    batch_size = 4

    np.random.seed(42)
    x_np = np.random.randn(batch_size, n_in).astype(np.float32)

    # JAX init and forward
    jax_block = JaxResnetBlock(n_out=n_out, n_hidden=n_hidden, training=False)
    key = jax.random.PRNGKey(0)
    jax_vars = jax_block.init(key, jnp.array(x_np))
    jax_out = np.array(jax_block.apply(jax_vars, jnp.array(x_np)))

    # Create PyTorch block and transfer weights
    torch_block = TorchResnetBlock(n_in=n_in, n_out=n_out, n_hidden=n_hidden)
    _transfer_resnet_block_weights(jax_vars["params"], torch_block)

    torch_block.eval()
    with torch.no_grad():
        torch_out = torch_block(torch.tensor(x_np)).numpy()

    max_diff = np.max(np.abs(jax_out - torch_out))
    print(f"ResnetBlock max diff: {max_diff:.2e}")
    np.testing.assert_allclose(jax_out, torch_out, atol=ATOL, rtol=RTOL)
    print("  PASSED")


def test_mlp_numerical():
    """Test MLP numerical equivalence with transferred weights."""
    from scvi.external.mrvi_jax._components import MLP as JaxMLP
    from scvi.external.mrvi_torch._components import MLP as TorchMLP

    n_in = 32
    n_out = 16
    n_hidden = 32
    batch_size = 4

    np.random.seed(42)
    x_np = np.random.randn(batch_size, n_in).astype(np.float32)

    jax_mlp = JaxMLP(n_out=n_out, n_hidden=n_hidden, n_layers=1, training=False)
    key = jax.random.PRNGKey(0)
    jax_vars = jax_mlp.init(key, jnp.array(x_np))
    jax_out = np.array(jax_mlp.apply(jax_vars, jnp.array(x_np)))

    torch_mlp = TorchMLP(n_in=n_in, n_out=n_out, n_hidden=n_hidden, n_layers=1)

    params = jax_vars["params"]
    _transfer_resnet_block_weights(params["ResnetBlock_0"], torch_mlp.resnet_blocks[0])
    _transfer_dense_weights(params["Dense_0"], torch_mlp.fc)

    torch_mlp.eval()
    with torch.no_grad():
        torch_out = torch_mlp(torch.tensor(x_np)).numpy()

    max_diff = np.max(np.abs(jax_out - torch_out))
    print(f"MLP max diff: {max_diff:.2e}")
    np.testing.assert_allclose(jax_out, torch_out, atol=ATOL, rtol=RTOL)
    print("  PASSED")


def test_conditional_normalization_numerical():
    """Test ConditionalNormalization numerical equivalence."""
    from scvi.external.mrvi_jax._components import (
        ConditionalNormalization as JaxCN,
    )
    from scvi.external.mrvi_torch._components import (
        ConditionalNormalization as TorchCN,
    )

    n_features = 32
    n_conditions = 5
    batch_size = 8

    np.random.seed(42)
    x_np = np.random.randn(batch_size, n_features).astype(np.float32)
    cond_np = np.random.randint(0, n_conditions, (batch_size, 1)).astype(np.float32)

    jax_cn = JaxCN(n_features=n_features, n_conditions=n_conditions, training=False)
    key = jax.random.PRNGKey(0)
    jax_vars = jax_cn.init(key, jnp.array(x_np), jnp.array(cond_np))
    jax_out = np.array(jax_cn.apply(jax_vars, jnp.array(x_np), jnp.array(cond_np)))

    torch_cn = TorchCN(n_features=n_features, n_conditions=n_conditions)
    params = jax_vars["params"]
    _transfer_embedding_weights(params["gamma_conditional"], torch_cn.gamma_embedding)
    _transfer_embedding_weights(params["beta_conditional"], torch_cn.beta_embedding)

    torch_cn.eval()
    with torch.no_grad():
        torch_out = torch_cn(torch.tensor(x_np), torch.tensor(cond_np), training=False).numpy()

    max_diff = np.max(np.abs(jax_out - torch_out))
    print(f"ConditionalNormalization max diff: {max_diff:.2e}")
    np.testing.assert_allclose(jax_out, torch_out, atol=ATOL, rtol=RTOL)
    print("  PASSED")


def test_attention_block_numerical():
    """Test AttentionBlock numerical equivalence with transferred weights."""
    from scvi.external.mrvi_jax._components import (
        AttentionBlock as JaxAttentionBlock,
    )
    from scvi.external.mrvi_torch._components import (
        AttentionBlock as TorchAttentionBlock,
    )

    query_dim = 10
    kv_dim = 16
    out_dim = 30
    outerprod_dim = 16
    n_channels = 4
    n_heads = 2
    batch_size = 4
    depth_per_head = n_channels  # qkv_features // n_heads

    np.random.seed(42)
    query_np = np.random.randn(batch_size, query_dim).astype(np.float32)
    kv_np = np.random.randn(batch_size, kv_dim).astype(np.float32)

    jax_block = JaxAttentionBlock(
        query_dim=query_dim,
        out_dim=out_dim,
        outerprod_dim=outerprod_dim,
        n_channels=n_channels,
        n_heads=n_heads,
        dropout_rate=0.0,
        n_hidden_mlp=32,
        n_layers_mlp=1,
        stop_gradients_mlp=False,
        training=False,
    )
    key = jax.random.PRNGKey(0)
    jax_vars = jax_block.init(
        {"params": key, "dropout": key},
        jnp.array(query_np),
        jnp.array(kv_np),
    )
    jax_out = np.array(
        jax_block.apply(
            jax_vars,
            jnp.array(query_np),
            jnp.array(kv_np),
            rngs={"dropout": key},
        )
    )

    torch_block = TorchAttentionBlock(
        query_dim=query_dim,
        kv_dim=kv_dim,
        out_dim=out_dim,
        outerprod_dim=outerprod_dim,
        n_channels=n_channels,
        n_heads=n_heads,
        dropout_rate=0.0,
        n_hidden_mlp=32,
        n_layers_mlp=1,
        stop_gradients_mlp=False,
    )

    params = jax_vars["params"]

    # Transfer query/kv projections (DenseGeneral → Linear)
    # JAX DenseGeneral((outerprod_dim, 1)) kernel: (query_dim, outerprod_dim, 1)
    # PyTorch Linear(query_dim, outerprod_dim) weight: (outerprod_dim, query_dim)
    query_kernel = np.array(params["DenseGeneral_0"]["kernel"])  # (query_dim, outerprod_dim, 1)
    torch_block.query_proj.weight.data = torch.tensor(query_kernel[:, :, 0].T)

    kv_kernel = np.array(params["DenseGeneral_1"]["kernel"])  # (kv_dim, outerprod_dim, 1)
    torch_block.kv_proj.weight.data = torch.tensor(kv_kernel[:, :, 0].T)

    # Transfer Q/K/V projections from MultiHeadDotProductAttention
    mha_params = params["MultiHeadDotProductAttention_0"]
    # Flax Q kernel: (1, n_heads, depth_per_head), bias: (n_heads, depth_per_head)
    # PyTorch Linear(1, n_heads*depth): weight (n_heads*depth, 1), bias (n_heads*depth)
    q_kernel = np.array(mha_params["query"]["kernel"])  # (1, n_heads, depth)
    q_bias = np.array(mha_params["query"]["bias"])  # (n_heads, depth)
    torch_block.q_proj.weight.data = torch.tensor(q_kernel.reshape(1, n_heads * depth_per_head).T)
    torch_block.q_proj.bias.data = torch.tensor(q_bias.reshape(-1))

    k_kernel = np.array(mha_params["key"]["kernel"])
    k_bias = np.array(mha_params["key"]["bias"])
    torch_block.k_proj.weight.data = torch.tensor(k_kernel.reshape(1, n_heads * depth_per_head).T)
    torch_block.k_proj.bias.data = torch.tensor(k_bias.reshape(-1))

    v_kernel = np.array(mha_params["value"]["kernel"])
    v_bias = np.array(mha_params["value"]["bias"])
    torch_block.v_proj.weight.data = torch.tensor(v_kernel.reshape(1, n_heads * depth_per_head).T)
    torch_block.v_proj.bias.data = torch.tensor(v_bias.reshape(-1))

    # Output projection
    # Flax: DenseGeneral(out_features, axis=(-2,-1)): kernel (n_heads, depth, out_features)
    # PyTorch: Linear(n_heads*depth, out_features): weight (out_features, n_heads*depth)
    out_kernel = np.array(mha_params["out"]["kernel"])  # (n_heads, depth_per_head, n_channels)
    out_bias = np.array(mha_params["out"]["bias"])  # (n_channels,)
    torch_block.out_proj.weight.data = torch.tensor(
        out_kernel.reshape(n_heads * depth_per_head, n_channels).T
    )
    torch_block.out_proj.bias.data = torch.tensor(out_bias)

    # Transfer MLP weights
    _transfer_mlp_weights(params["MLP_0"], torch_block.mlp_eps)
    _transfer_mlp_weights(params["MLP_1"], torch_block.mlp_residual)

    torch_block.eval()
    with torch.no_grad():
        torch_out = torch_block(torch.tensor(query_np), torch.tensor(kv_np)).numpy()

    max_diff = np.max(np.abs(jax_out - torch_out))
    print(f"AttentionBlock max diff: {max_diff:.2e}")
    # Use compound tolerance because softmax numerical differences (~2e-7)
    # amplify through downstream MLPs in float32
    np.testing.assert_allclose(jax_out, torch_out, atol=ATOL_COMPOUND, rtol=RTOL_COMPOUND)
    print("  PASSED")


def test_normal_dist_output_nn_numerical():
    """Test NormalDistOutputNN numerical equivalence."""
    from scvi.external.mrvi_jax._components import (
        NormalDistOutputNN as JaxNDONN,
    )
    from scvi.external.mrvi_torch._components import (
        NormalDistOutputNN as TorchNDONN,
    )

    n_in = 32
    n_out = 10
    n_hidden = 32
    batch_size = 4

    np.random.seed(42)
    x_np = np.random.randn(batch_size, n_in).astype(np.float32)

    jax_nn = JaxNDONN(n_out=n_out, n_hidden=n_hidden, n_layers=1, training=False)
    key = jax.random.PRNGKey(0)
    jax_vars = jax_nn.init(key, jnp.array(x_np))
    jax_dist = jax_nn.apply(jax_vars, jnp.array(x_np))
    jax_mean = np.array(jax_dist.mean)
    jax_scale = np.array(jax_dist.scale)

    torch_nn = TorchNDONN(n_in=n_in, n_out=n_out, n_hidden=n_hidden, n_layers=1)

    params = jax_vars["params"]
    _transfer_resnet_block_weights(params["ResnetBlock_0"], torch_nn.resnet_blocks[0])
    _transfer_dense_weights(params["Dense_0"], torch_nn.fc_mean)
    _transfer_dense_weights(params["Dense_1"], torch_nn.fc_scale[0])

    torch_nn.eval()
    with torch.no_grad():
        torch_dist = torch_nn(torch.tensor(x_np))
        torch_mean = torch_dist.mean.numpy()
        torch_scale = torch_dist.scale.numpy()

    mean_diff = np.max(np.abs(jax_mean - torch_mean))
    scale_diff = np.max(np.abs(jax_scale - torch_scale))
    print(f"NormalDistOutputNN mean max diff: {mean_diff:.2e}, scale max diff: {scale_diff:.2e}")
    np.testing.assert_allclose(jax_mean, torch_mean, atol=ATOL, rtol=RTOL)
    np.testing.assert_allclose(jax_scale, torch_scale, atol=ATOL, rtol=RTOL)
    print("  PASSED")


if __name__ == "__main__":
    print("=" * 60)
    print("NUMERICAL EQUIVALENCE TESTS")
    print("=" * 60)

    print("\n--- ResnetBlock ---")
    test_resnet_block_numerical()

    print("\n--- MLP ---")
    test_mlp_numerical()

    print("\n--- ConditionalNormalization ---")
    test_conditional_normalization_numerical()

    print("\n--- NormalDistOutputNN ---")
    test_normal_dist_output_nn_numerical()

    print("\n--- AttentionBlock ---")
    test_attention_block_numerical()

    print("\n" + "=" * 60)
    print("ALL NUMERICAL EQUIVALENCE TESTS PASSED!")
    print("=" * 60)
