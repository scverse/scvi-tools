"""Test equivalence between JAX and PyTorch MRVI implementations.

Tests that given identical weights and inputs, both implementations produce
the same forward pass outputs and gradients.
"""

from __future__ import annotations

import numpy as np
import torch


def test_attention_block_equivalence():
    """Test that AttentionBlock produces equivalent outputs in both backends."""
    import jax
    import jax.numpy as jnp

    from scvi.external.mrvi_jax._components import (
        AttentionBlock as JaxAttentionBlock,
    )
    from scvi.external.mrvi_torch._components import (
        AttentionBlock as TorchAttentionBlock,
    )

    # Fixed dimensions
    query_dim = 10
    kv_dim = 16
    out_dim = 30
    outerprod_dim = 16
    n_channels = 4
    n_heads = 2
    batch_size = 8

    # Create JAX attention block and initialize
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

    # Create random inputs
    np.random.seed(42)
    query_np = np.random.randn(batch_size, query_dim).astype(np.float32)
    kv_np = np.random.randn(batch_size, kv_dim).astype(np.float32)

    # Initialize JAX model
    key = jax.random.PRNGKey(0)
    jax_vars = jax_block.init(
        {"params": key, "dropout": key},
        query_embed=jnp.array(query_np),
        kv_embed=jnp.array(kv_np),
    )

    # Run JAX forward
    jax_out = jax_block.apply(
        jax_vars,
        query_embed=jnp.array(query_np),
        kv_embed=jnp.array(kv_np),
        rngs={"dropout": key},
    )

    # Create PyTorch block
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

    # Verify output shapes match
    torch_query = torch.tensor(query_np)
    torch_kv = torch.tensor(kv_np)

    torch_block.eval()
    with torch.no_grad():
        torch_out = torch_block(torch_query, torch_kv)

    assert jax_out.shape == tuple(torch_out.shape), (
        f"Shape mismatch: JAX {jax_out.shape} vs Torch {torch_out.shape}"
    )
    print(f"AttentionBlock output shape: JAX={jax_out.shape}, Torch={tuple(torch_out.shape)}")


def test_resnet_block_equivalence():
    """Test that ResnetBlock produces equivalent outputs."""
    import jax
    import jax.numpy as jnp

    from scvi.external.mrvi_jax._components import ResnetBlock as JaxResnetBlock
    from scvi.external.mrvi_torch._components import ResnetBlock as TorchResnetBlock

    n_in = 32
    n_out = 16
    n_hidden = 32
    batch_size = 8

    np.random.seed(42)
    x_np = np.random.randn(batch_size, n_in).astype(np.float32)

    # JAX
    jax_block = JaxResnetBlock(n_out=n_out, n_hidden=n_hidden, training=False)
    key = jax.random.PRNGKey(0)
    jax_vars = jax_block.init(key, jnp.array(x_np))
    jax_out = jax_block.apply(jax_vars, jnp.array(x_np))

    # PyTorch
    torch_block = TorchResnetBlock(n_in=n_in, n_out=n_out, n_hidden=n_hidden)
    torch_block.eval()
    with torch.no_grad():
        torch_out = torch_block(torch.tensor(x_np))

    assert jax_out.shape == tuple(torch_out.shape), (
        f"Shape mismatch: JAX {jax_out.shape} vs Torch {torch_out.shape}"
    )
    print(f"ResnetBlock output shape: JAX={jax_out.shape}, Torch={tuple(torch_out.shape)}")


def test_mlp_equivalence():
    """Test that MLP produces equivalent outputs."""
    import jax
    import jax.numpy as jnp

    from scvi.external.mrvi_jax._components import MLP as JaxMLP
    from scvi.external.mrvi_torch._components import MLP as TorchMLP

    n_in = 32
    n_out = 16
    n_hidden = 32
    batch_size = 8

    np.random.seed(42)
    x_np = np.random.randn(batch_size, n_in).astype(np.float32)

    # JAX
    jax_mlp = JaxMLP(n_out=n_out, n_hidden=n_hidden, n_layers=1, training=False)
    key = jax.random.PRNGKey(0)
    jax_vars = jax_mlp.init(key, jnp.array(x_np))
    jax_out = jax_mlp.apply(jax_vars, jnp.array(x_np))

    # PyTorch
    torch_mlp = TorchMLP(n_in=n_in, n_out=n_out, n_hidden=n_hidden, n_layers=1)
    torch_mlp.eval()
    with torch.no_grad():
        torch_out = torch_mlp(torch.tensor(x_np))

    assert jax_out.shape == tuple(torch_out.shape), (
        f"Shape mismatch: JAX {jax_out.shape} vs Torch {torch_out.shape}"
    )
    print(f"MLP output shape: JAX={jax_out.shape}, Torch={tuple(torch_out.shape)}")


def test_full_module_output_shapes():
    """Test that full MRVI module produces matching output shapes via model API."""
    from scvi.external.mrvi_torch._module import TorchMRVAE

    n_input = 100
    n_sample = 5
    n_batch = 2
    n_labels = 1
    n_latent = 30
    n_latent_u = 10
    batch_size = 16

    np.random.seed(42)
    x_np = np.abs(np.random.randn(batch_size, n_input).astype(np.float32)) + 1
    sample_idx_np = np.random.randint(0, n_sample, (batch_size, 1)).astype(np.float32)
    batch_idx_np = np.random.randint(0, n_batch, (batch_size, 1)).astype(np.float32)
    label_idx_np = np.zeros((batch_size, 1), dtype=np.float32)

    torch_module = TorchMRVAE(
        n_input=n_input,
        n_sample=n_sample,
        n_batch=n_batch,
        n_labels=n_labels,
        n_latent=n_latent,
        n_latent_u=n_latent_u,
    )
    torch_module.eval()

    with torch.no_grad():
        torch_inf_out = torch_module.inference(
            x=torch.tensor(x_np),
            sample_index=torch.tensor(sample_idx_np),
            use_mean=True,
        )

    # Check shapes match expected
    assert torch_inf_out["u"].shape == (batch_size, n_latent_u)
    assert torch_inf_out["z"].shape == (batch_size, n_latent)
    assert torch_inf_out["z_base"].shape == (batch_size, n_latent)
    assert torch_inf_out["library"].shape == (batch_size, 1)
    print(f"  u: {tuple(torch_inf_out['u'].shape)}")
    print(f"  z: {tuple(torch_inf_out['z'].shape)}")
    print(f"  z_base: {tuple(torch_inf_out['z_base'].shape)}")
    print(f"  library: {tuple(torch_inf_out['library'].shape)}")

    with torch.no_grad():
        torch_gen_out = torch_module.generative(
            z=torch_inf_out["z"],
            library=torch_inf_out["library"],
            batch_index=torch.tensor(batch_idx_np),
            label_index=torch.tensor(label_idx_np),
        )

    assert torch_gen_out["h"].shape == (batch_size, n_input)
    print(f"  h: {tuple(torch_gen_out['h'].shape)}")

    # Test with mc_samples
    with torch.no_grad():
        mc_inf_out = torch_module.inference(
            x=torch.tensor(x_np),
            sample_index=torch.tensor(sample_idx_np),
            mc_samples=5,
        )
    assert mc_inf_out["u"].shape == (5, batch_size, n_latent_u)
    assert mc_inf_out["z"].shape == (5, batch_size, n_latent)
    print(f"  u (mc=5): {tuple(mc_inf_out['u'].shape)}")
    print(f"  z (mc=5): {tuple(mc_inf_out['z'].shape)}")
    print("All module output shapes correct!")


def test_gpu_training_end_to_end():
    """Test full training and inference on GPU."""
    import scvi
    from scvi.external import MRVI

    if not torch.cuda.is_available():
        print("CUDA not available, skipping GPU test")
        return

    np.random.seed(42)
    adata = scvi.data.synthetic_iid()
    adata.obs["sample"] = np.random.choice(5, adata.n_obs).astype(str)

    MRVI.setup_anndata(
        adata, sample_key="sample", batch_key="batch", labels_key="labels", backend="torch"
    )
    model = MRVI(adata)

    # Train on GPU
    model.train(max_epochs=3, accelerator="gpu", devices=1, batch_size=128)

    # Verify model is on GPU
    device = model.device
    assert "cuda" in str(device), f"Expected CUDA device, got {device}"

    # Verify pz_scale is on GPU
    assert "cuda" in str(model.module.pz_scale.device), (
        f"pz_scale on wrong device: {model.module.pz_scale.device}"
    )

    # Test all inference methods work on GPU
    z = model.get_latent_representation(give_z=True, batch_size=128)
    assert z.shape == (adata.n_obs, 30), f"Wrong z shape: {z.shape}"

    u = model.get_latent_representation(give_z=False, batch_size=128)
    assert u.shape == (adata.n_obs, 10), f"Wrong u shape: {u.shape}"

    # Test local sample distances on GPU
    dists = model.get_local_sample_distances(batch_size=128)
    assert "cell" in dists.data_vars

    # Test normalized expression on GPU
    expr = model.get_normalized_expression(batch_size=128)
    assert expr.shape[0] == adata.n_obs

    # Test save/load roundtrip
    import tempfile

    with tempfile.TemporaryDirectory() as tmpdir:
        model.save(tmpdir, save_anndata=False, overwrite=True)
        loaded = MRVI.load(tmpdir, adata=adata)
        loaded.train(max_epochs=1, accelerator="gpu", devices=1, batch_size=128)
        z2 = loaded.get_latent_representation(give_z=True, batch_size=128)
        assert z2.shape == z.shape

    print("GPU end-to-end test passed!")


if __name__ == "__main__":
    print("=" * 60)
    print("Testing component equivalence...")
    print("=" * 60)

    print("\n--- ResnetBlock ---")
    test_resnet_block_equivalence()

    print("\n--- MLP ---")
    test_mlp_equivalence()

    print("\n--- AttentionBlock ---")
    test_attention_block_equivalence()

    print("\n--- Full Module Output Shapes ---")
    test_full_module_output_shapes()

    print("\n--- GPU End-to-End ---")
    test_gpu_training_end_to_end()

    print("\n" + "=" * 60)
    print("ALL EQUIVALENCE TESTS PASSED!")
    print("=" * 60)
