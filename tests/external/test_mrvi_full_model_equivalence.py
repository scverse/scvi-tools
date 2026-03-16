"""Full model numerical equivalence test between JAX and PyTorch MRVI.

Initializes the JAX module, transfers every parameter to PyTorch,
runs both with identical inputs using use_mean=True (deterministic),
and compares every intermediate output: encoder XU, encoder UZ,
decoder, generative, and loss.
"""

from __future__ import annotations

import os

os.environ["JAX_DEFAULT_MATMUL_PRECISION"] = "float32"
os.environ["JAX_PLATFORMS"] = "cpu"

import numpy as np
import torch

torch.backends.cuda.matmul.allow_tf32 = False
torch.backends.cudnn.allow_tf32 = False
torch.set_float32_matmul_precision("highest")

import jax  # noqa: E402
import jax.numpy as jnp  # noqa: E402

jax.config.update("jax_default_matmul_precision", "float32")

from scvi import REGISTRY_KEYS  # noqa: E402
from scvi.external.mrvi_jax._module import JaxMRVAE  # noqa: E402
from scvi.external.mrvi_torch._module import TorchMRVAE  # noqa: E402

# Tolerances: simple paths stay tight, paths through attention are looser
ATOL_SIMPLE = 1e-5  # e.g. encoder XU outputs (no attention)
ATOL_ATTN = 5e-4  # paths that go through attention (softmax amplification)
RTOL = 1e-3


# ── weight transfer helpers ──────────────────────────────────────────────────


def _t(arr):
    """Convert JAX/numpy array to float32 torch tensor."""
    return torch.tensor(np.array(arr)).float()


def _dense(jax_p, torch_linear):
    """Transfer Dense/DenseGeneral → nn.Linear."""
    torch_linear.weight.data = _t(jax_p["kernel"].T)
    if "bias" in jax_p:
        torch_linear.bias.data = _t(jax_p["bias"])


def _ln(jax_p, torch_ln):
    """Transfer LayerNorm."""
    if "scale" in jax_p:
        torch_ln.weight.data = _t(jax_p["scale"])
    if "bias" in jax_p:
        torch_ln.bias.data = _t(jax_p["bias"])


def _embed(jax_p, torch_emb):
    """Transfer Embedding."""
    torch_emb.weight.data = _t(jax_p["embedding"])


def _resnet_block(jp, tb):
    """Transfer ResnetBlock."""
    _dense(jp["Dense_0"], tb.fc1)
    _ln(jp["LayerNorm_0"], tb.layer_norm1)
    if tb.fc_match is not None:
        _dense(jp["Dense_1"], tb.fc_match)
        _dense(jp["Dense_2"], tb.fc2)
    else:
        _dense(jp["Dense_1"], tb.fc2)
    _ln(jp["LayerNorm_1"], tb.layer_norm2)


def _mlp(jp, tm):
    """Transfer MLP (ResnetBlocks + final Dense)."""
    for i in range(len(tm.resnet_blocks)):
        _resnet_block(jp[f"ResnetBlock_{i}"], tm.resnet_blocks[i])
    _dense(jp["Dense_0"], tm.fc)


def _normal_dist_nn(jp, tn):
    """Transfer NormalDistOutputNN."""
    for i, rb in enumerate(tn.resnet_blocks):
        _resnet_block(jp[f"ResnetBlock_{i}"], rb)
    _dense(jp["Dense_0"], tn.fc_mean)
    _dense(jp["Dense_1"], tn.fc_scale[0])


def _attention_block(jp, tab):
    """Transfer AttentionBlock (projections + MHDPA + MLPs)."""
    n_heads = tab.n_heads
    depth = tab.depth_per_head
    n_channels = tab.n_channels

    # DenseGeneral query/kv projections → Linear
    tab.query_proj.weight.data = _t(jp["DenseGeneral_0"]["kernel"][:, :, 0].T)
    tab.kv_proj.weight.data = _t(jp["DenseGeneral_1"]["kernel"][:, :, 0].T)

    # MHDPA Q/K/V projections
    mha = jp["MultiHeadDotProductAttention_0"]
    for name, proj in [("query", tab.q_proj), ("key", tab.k_proj), ("value", tab.v_proj)]:
        k = np.array(mha[name]["kernel"])  # (1, n_heads, depth)
        b = np.array(mha[name]["bias"])  # (n_heads, depth)
        proj.weight.data = _t(k.reshape(1, n_heads * depth).T)
        proj.bias.data = _t(b.reshape(-1))

    # MHDPA output projection
    ok = np.array(mha["out"]["kernel"])  # (n_heads, depth, n_channels)
    ob = np.array(mha["out"]["bias"])  # (n_channels,)
    tab.out_proj.weight.data = _t(ok.reshape(n_heads * depth, n_channels).T)
    tab.out_proj.bias.data = _t(ob)

    # MLPs
    _mlp(jp["MLP_0"], tab.mlp_eps)
    _mlp(jp["MLP_1"], tab.mlp_residual)


# ── full module transfer ─────────────────────────────────────────────────────


def transfer_all_weights(jax_params, torch_module):
    """Transfer every parameter from JAX MRVI module to PyTorch MRVI module."""
    p = jax_params

    # ── EncoderXU (qu) ──
    tqu = torch_module.qu
    _dense(p["qu"]["Dense_0"], tqu.fc1)
    _dense(p["qu"]["Dense_1"], tqu.fc2)
    _embed(
        p["qu"]["ConditionalNormalization_0"]["gamma_conditional"],
        tqu.conditional_norm1.gamma_embedding,
    )
    _embed(
        p["qu"]["ConditionalNormalization_0"]["beta_conditional"],
        tqu.conditional_norm1.beta_embedding,
    )
    _embed(
        p["qu"]["ConditionalNormalization_1"]["gamma_conditional"],
        tqu.conditional_norm2.gamma_embedding,
    )
    _embed(
        p["qu"]["ConditionalNormalization_1"]["beta_conditional"],
        tqu.conditional_norm2.beta_embedding,
    )
    _embed(p["qu"]["Embed_0"], tqu.sample_embed)
    _normal_dist_nn(p["qu"]["NormalDistOutputNN_0"], tqu.output_nn)

    # ── EncoderUZ (qz) ──
    tqz = torch_module.qz
    _ln(p["qz"]["u_ln"], tqz.layer_norm)
    _embed(p["qz"]["Embed_0"], tqz.embedding)
    _ln(p["qz"]["sample_embed_ln"], tqz.layer_norm_embed)
    _attention_block(p["qz"]["AttentionBlock_0"], tqz.attention_block)
    # fc projection (u → z_base)
    if tqz.fc is not None:
        _dense(p["qz"]["Dense_0"], tqz.fc)

    # ── DecoderZXAttention (px) ──
    tpx = torch_module.px
    _ln(p["px"]["u_ln"], tpx.layer_norm)
    _embed(p["px"]["Embed_0"], tpx.batch_embedding)
    _ln(p["px"]["batch_embed_ln"], tpx.layer_norm_batch_embed)
    _attention_block(p["px"]["AttentionBlock_0"], tpx.attention_block)
    _dense(p["px"]["Dense_0"], tpx.fc)
    tpx.px_r.data = _t(p["px"]["px_r"])

    # ── Mixture prior parameters ──
    if torch_module.u_prior_mixture:
        torch_module.u_prior_logits.data = _t(p["u_prior_logits"])
        torch_module.u_prior_means.data = _t(p["u_prior_means"])
        torch_module.u_prior_scales.data = _t(p["u_prior_scales"])

    # ── pz_scale ──
    if hasattr(torch_module, "pz_scale"):
        if isinstance(torch_module.pz_scale, torch.nn.Parameter):
            torch_module.pz_scale.data = _t(p["pz_scale"])
        else:
            # buffer
            pz_jax = (
                np.array(p["pz_scale"]) if "pz_scale" in p else np.zeros(torch_module.n_latent)
            )
            torch_module.pz_scale.copy_(_t(pz_jax))


# ── comparison helper ────────────────────────────────────────────────────────


def _compare(name, jax_val, torch_val, atol=ATOL_SIMPLE):
    """Compare and print result."""
    j = np.array(jax_val)
    t = (
        torch_val.detach().cpu().numpy()
        if isinstance(torch_val, torch.Tensor)
        else np.array(torch_val)
    )
    diff = np.max(np.abs(j - t))
    status = "OK" if diff < atol else "FAIL"
    print(f"  {name:40s}  shape {str(j.shape):20s}  max_diff={diff:.2e}  [{status}]")
    assert diff < atol, f"{name}: max_diff={diff:.2e} exceeds atol={atol:.0e}"
    return diff


# ── main test ────────────────────────────────────────────────────────────────


def test_full_model_equivalence():
    """End-to-end: init JAX → transfer weights → compare every intermediate."""
    n_input, n_sample, n_batch, n_labels = 100, 5, 2, 1
    n_latent, n_latent_u = 30, 10
    bs = 16

    # ── create deterministic inputs ──
    np.random.seed(0)
    x_np = (np.abs(np.random.randn(bs, n_input)) + 0.1).astype(np.float32)
    sample_np = np.random.randint(0, n_sample, (bs, 1)).astype(np.float32)
    batch_np = np.random.randint(0, n_batch, (bs, 1)).astype(np.float32)
    label_np = np.zeros((bs, 1), dtype=np.float32)

    # ── JAX module init ──
    jax_mod = JaxMRVAE(
        n_input=n_input,
        n_sample=n_sample,
        n_batch=n_batch,
        n_labels=n_labels,
        n_latent=n_latent,
        n_latent_u=n_latent_u,
        training=False,
    )

    tensors_jax = {
        REGISTRY_KEYS.X_KEY: jnp.array(x_np),
        REGISTRY_KEYS.SAMPLE_KEY: jnp.array(sample_np),
        REGISTRY_KEYS.BATCH_KEY: jnp.array(batch_np),
        REGISTRY_KEYS.LABELS_KEY: jnp.array(label_np),
    }

    key = jax.random.PRNGKey(42)
    keys = jax.random.split(key, 4)
    rngs = {"params": keys[0], "u": keys[1], "dropout": keys[2], "eps": keys[3]}

    variables = jax_mod.init(rngs, tensors_jax)
    jax_params = variables["params"]

    # ── PyTorch module init + weight transfer ──
    torch_mod = TorchMRVAE(
        n_input=n_input,
        n_sample=n_sample,
        n_batch=n_batch,
        n_labels=n_labels,
        n_latent=n_latent,
        n_latent_u=n_latent_u,
    )
    torch_mod.eval()
    transfer_all_weights(jax_params, torch_mod)

    # ── 1. Inference (use_mean=True → deterministic, no sampling) ──
    print("\n=== INFERENCE (use_mean=True) ===")

    jax_inf = jax_mod.apply(
        variables,
        x=jnp.array(x_np),
        sample_index=jnp.array(sample_np),
        use_mean=True,
        method=jax_mod.inference,
        rngs=rngs,
    )

    with torch.no_grad():
        torch_inf = torch_mod.inference(
            x=torch.tensor(x_np),
            sample_index=torch.tensor(sample_np),
            use_mean=True,
        )

    _compare("qu.mean", jax_inf["qu"].mean, torch_inf["qu"].mean, ATOL_SIMPLE)
    _compare("qu.scale", jax_inf["qu"].scale, torch_inf["qu"].scale, ATOL_SIMPLE)
    _compare("u", jax_inf["u"], torch_inf["u"], ATOL_SIMPLE)
    _compare("z_base", jax_inf["z_base"], torch_inf["z_base"], ATOL_ATTN)
    _compare("eps (attention residual)", jax_inf["eps"], torch_inf["eps"], ATOL_ATTN)
    _compare("z = z_base + eps", jax_inf["z"], torch_inf["z"], ATOL_ATTN)
    _compare("library", jax_inf["library"], torch_inf["library"], ATOL_SIMPLE)

    # ── 2. Generative ──
    print("\n=== GENERATIVE ===")

    jax_gen = jax_mod.apply(
        variables,
        z=jax_inf["z"],
        library=jax_inf["library"],
        batch_index=jnp.array(batch_np),
        label_index=jnp.array(label_np),
        method=jax_mod.generative,
        rngs=rngs,
    )

    with torch.no_grad():
        torch_gen = torch_mod.generative(
            z=torch_inf["z"],
            library=torch_inf["library"],
            batch_index=torch.tensor(batch_np),
            label_index=torch.tensor(label_np),
        )

    _compare("h (normalized expr)", jax_gen["h"], torch_gen["h"], ATOL_ATTN)
    _compare("px.mean", jax_gen["px"].mean, torch_gen["px"].mean, ATOL_ATTN)

    # ── 3. Loss ──
    print("\n=== LOSS ===")

    tensors_torch = {
        REGISTRY_KEYS.X_KEY: torch.tensor(x_np),
        REGISTRY_KEYS.SAMPLE_KEY: torch.tensor(sample_np),
        REGISTRY_KEYS.BATCH_KEY: torch.tensor(batch_np),
        REGISTRY_KEYS.LABELS_KEY: torch.tensor(label_np),
    }

    jax_loss_out = jax_mod.apply(
        variables,
        tensors=tensors_jax,
        inference_outputs=jax_inf,
        generative_outputs=jax_gen,
        kl_weight=1.0,
        method=jax_mod.loss,
        rngs=rngs,
    )

    with torch.no_grad():
        torch_loss_out = torch_mod.loss(
            tensors=tensors_torch,
            inference_outputs=torch_inf,
            generative_outputs=torch_gen,
            kl_weight=1.0,
        )

    # Both LossOutputs may wrap per-obs losses in dicts; extract the raw tensors
    def _extract_loss(val):
        if isinstance(val, dict):
            return sum(val.values())
        return val

    _compare(
        "reconstruction_loss",
        _extract_loss(jax_loss_out.reconstruction_loss),
        _extract_loss(torch_loss_out.reconstruction_loss),
        ATOL_ATTN,
    )
    _compare(
        "kl_local",
        _extract_loss(jax_loss_out.kl_local),
        _extract_loss(torch_loss_out.kl_local),
        ATOL_ATTN,
    )
    _compare("total loss", jax_loss_out.loss, torch_loss_out.loss, ATOL_ATTN)

    print("\n" + "=" * 72)
    print("FULL MODEL EQUIVALENCE TEST PASSED")
    print("=" * 72)


# ── gradient comparison helpers ──────────────────────────────────────────────


def _compare_grad(name, jax_grad_val, torch_param, atol):
    """Compare a JAX gradient array against a PyTorch param's .grad, with layout transform."""
    j = np.array(jax_grad_val)
    assert torch_param.grad is not None, f"{name}: PyTorch grad is None"
    t = torch_param.grad.detach().cpu().numpy()
    assert j.shape == t.shape, f"{name}: shape mismatch JAX {j.shape} vs Torch {t.shape}"
    diff = np.max(np.abs(j - t))
    status = "OK" if diff < atol else "FAIL"
    print(f"  grad {name:40s}  shape {str(j.shape):20s}  max_diff={diff:.2e}  [{status}]")
    assert diff < atol, f"grad {name}: max_diff={diff:.2e} exceeds atol={atol:.0e}"


def _compare_dense_grad(name, jax_grad, torch_linear, atol):
    """Compare Dense/Linear gradients (kernel needs transpose)."""
    _compare_grad(f"{name}.weight", np.array(jax_grad["kernel"]).T, torch_linear.weight, atol)
    if "bias" in jax_grad and torch_linear.bias is not None:
        _compare_grad(f"{name}.bias", jax_grad["bias"], torch_linear.bias, atol)


def _compare_ln_grad(name, jax_grad, torch_ln, atol):
    """Compare LayerNorm gradients."""
    if "scale" in jax_grad:
        _compare_grad(f"{name}.weight", jax_grad["scale"], torch_ln.weight, atol)
    if "bias" in jax_grad:
        _compare_grad(f"{name}.bias", jax_grad["bias"], torch_ln.bias, atol)


def _compare_embed_grad(name, jax_grad, torch_emb, atol):
    """Compare Embedding gradients."""
    _compare_grad(f"{name}.weight", jax_grad["embedding"], torch_emb.weight, atol)


def _compare_resnet_block_grad(name, jg, tb, atol):
    """Compare ResnetBlock gradients."""
    _compare_dense_grad(f"{name}.fc1", jg["Dense_0"], tb.fc1, atol)
    _compare_ln_grad(f"{name}.ln1", jg["LayerNorm_0"], tb.layer_norm1, atol)
    if tb.fc_match is not None:
        _compare_dense_grad(f"{name}.fc_match", jg["Dense_1"], tb.fc_match, atol)
        _compare_dense_grad(f"{name}.fc2", jg["Dense_2"], tb.fc2, atol)
    else:
        _compare_dense_grad(f"{name}.fc2", jg["Dense_1"], tb.fc2, atol)
    _compare_ln_grad(f"{name}.ln2", jg["LayerNorm_1"], tb.layer_norm2, atol)


def _compare_mlp_grad(name, jg, tm, atol):
    """Compare MLP gradients."""
    for i in range(len(tm.resnet_blocks)):
        _compare_resnet_block_grad(
            f"{name}.rb{i}", jg[f"ResnetBlock_{i}"], tm.resnet_blocks[i], atol
        )
    _compare_dense_grad(f"{name}.fc", jg["Dense_0"], tm.fc, atol)


def _compare_normal_dist_nn_grad(name, jg, tn, atol):
    """Compare NormalDistOutputNN gradients."""
    for i, rb in enumerate(tn.resnet_blocks):
        _compare_resnet_block_grad(f"{name}.rb{i}", jg[f"ResnetBlock_{i}"], rb, atol)
    _compare_dense_grad(f"{name}.fc_mean", jg["Dense_0"], tn.fc_mean, atol)
    _compare_dense_grad(f"{name}.fc_scale", jg["Dense_1"], tn.fc_scale[0], atol)


def _compare_attention_block_grad(name, jg, tab, atol):
    """Compare AttentionBlock gradients (with proper reshape for MHDPA weights)."""
    n_heads = tab.n_heads
    depth = tab.depth_per_head
    n_channels = tab.n_channels

    # DenseGeneral projections: kernel → weight (transposed)
    _compare_grad(
        f"{name}.query_proj.w",
        np.array(jg["DenseGeneral_0"]["kernel"])[:, :, 0].T,
        tab.query_proj.weight,
        atol,
    )
    _compare_grad(
        f"{name}.kv_proj.w",
        np.array(jg["DenseGeneral_1"]["kernel"])[:, :, 0].T,
        tab.kv_proj.weight,
        atol,
    )

    # MHDPA Q/K/V: kernel (1, n_heads, depth) → weight (n_heads*depth, 1)
    mha_g = jg["MultiHeadDotProductAttention_0"]
    for qkv_name, proj in [("query", tab.q_proj), ("key", tab.k_proj), ("value", tab.v_proj)]:
        kg = np.array(mha_g[qkv_name]["kernel"])
        bg = np.array(mha_g[qkv_name]["bias"])
        _compare_grad(f"{name}.{qkv_name}.w", kg.reshape(1, n_heads * depth).T, proj.weight, atol)
        _compare_grad(f"{name}.{qkv_name}.b", bg.reshape(-1), proj.bias, atol)

    # Output projection: kernel (n_heads, depth, n_channels) → weight (n_channels, n_heads*depth)
    okg = np.array(mha_g["out"]["kernel"])
    obg = np.array(mha_g["out"]["bias"])
    _compare_grad(
        f"{name}.out.w", okg.reshape(n_heads * depth, n_channels).T, tab.out_proj.weight, atol
    )
    _compare_grad(f"{name}.out.b", obg, tab.out_proj.bias, atol)

    # MLPs
    _compare_mlp_grad(f"{name}.mlp_eps", jg["MLP_0"], tab.mlp_eps, atol)
    _compare_mlp_grad(f"{name}.mlp_res", jg["MLP_1"], tab.mlp_residual, atol)


def compare_all_gradients(jax_grads, torch_mod, atol_simple, atol_attn):
    """Compare every gradient from JAX grad tree against PyTorch .grad values."""
    g = jax_grads

    print("\n  -- EncoderXU (qu) --")
    tqu = torch_mod.qu
    _compare_dense_grad("qu.fc1", g["qu"]["Dense_0"], tqu.fc1, atol_simple)
    _compare_dense_grad("qu.fc2", g["qu"]["Dense_1"], tqu.fc2, atol_simple)
    _compare_embed_grad(
        "qu.cn0.gamma",
        g["qu"]["ConditionalNormalization_0"]["gamma_conditional"],
        tqu.conditional_norm1.gamma_embedding,
        atol_simple,
    )
    _compare_embed_grad(
        "qu.cn0.beta",
        g["qu"]["ConditionalNormalization_0"]["beta_conditional"],
        tqu.conditional_norm1.beta_embedding,
        atol_simple,
    )
    _compare_embed_grad(
        "qu.cn1.gamma",
        g["qu"]["ConditionalNormalization_1"]["gamma_conditional"],
        tqu.conditional_norm2.gamma_embedding,
        atol_simple,
    )
    _compare_embed_grad(
        "qu.cn1.beta",
        g["qu"]["ConditionalNormalization_1"]["beta_conditional"],
        tqu.conditional_norm2.beta_embedding,
        atol_simple,
    )
    _compare_embed_grad("qu.sample_embed", g["qu"]["Embed_0"], tqu.sample_embed, atol_simple)
    _compare_normal_dist_nn_grad(
        "qu.output_nn", g["qu"]["NormalDistOutputNN_0"], tqu.output_nn, atol_simple
    )

    print("\n  -- EncoderUZ (qz) --")
    tqz = torch_mod.qz
    _compare_ln_grad("qz.u_ln", g["qz"]["u_ln"], tqz.layer_norm, atol_attn)
    _compare_embed_grad("qz.embed", g["qz"]["Embed_0"], tqz.embedding, atol_attn)
    _compare_ln_grad("qz.embed_ln", g["qz"]["sample_embed_ln"], tqz.layer_norm_embed, atol_attn)
    _compare_attention_block_grad(
        "qz.attn", g["qz"]["AttentionBlock_0"], tqz.attention_block, atol_attn
    )
    if tqz.fc is not None:
        _compare_dense_grad("qz.fc", g["qz"]["Dense_0"], tqz.fc, atol_attn)

    print("\n  -- DecoderZXAttention (px) --")
    tpx = torch_mod.px
    _compare_ln_grad("px.u_ln", g["px"]["u_ln"], tpx.layer_norm, atol_attn)
    _compare_embed_grad("px.batch_embed", g["px"]["Embed_0"], tpx.batch_embedding, atol_attn)
    _compare_ln_grad(
        "px.batch_ln", g["px"]["batch_embed_ln"], tpx.layer_norm_batch_embed, atol_attn
    )
    _compare_attention_block_grad(
        "px.attn", g["px"]["AttentionBlock_0"], tpx.attention_block, atol_attn
    )
    _compare_dense_grad("px.fc", g["px"]["Dense_0"], tpx.fc, atol_attn)
    _compare_grad("px.px_r", g["px"]["px_r"], tpx.px_r, atol_attn)

    print("\n  -- Prior parameters --")
    if torch_mod.u_prior_mixture:
        _compare_grad("u_prior_logits", g["u_prior_logits"], torch_mod.u_prior_logits, atol_simple)
        _compare_grad("u_prior_means", g["u_prior_means"], torch_mod.u_prior_means, atol_simple)
        _compare_grad("u_prior_scales", g["u_prior_scales"], torch_mod.u_prior_scales, atol_simple)


# ── gradient test ────────────────────────────────────────────────────────────


def test_gradient_equivalence():
    """Compare gradients from a single forward-backward pass."""
    n_input, n_sample, n_batch, n_labels = 100, 5, 2, 1
    n_latent, n_latent_u = 30, 10
    bs = 16

    np.random.seed(0)
    x_np = (np.abs(np.random.randn(bs, n_input)) + 0.1).astype(np.float32)
    sample_np = np.random.randint(0, n_sample, (bs, 1)).astype(np.float32)
    batch_np = np.random.randint(0, n_batch, (bs, 1)).astype(np.float32)
    label_np = np.zeros((bs, 1), dtype=np.float32)

    # ── JAX: init, transfer, compute gradients ──
    jax_mod = JaxMRVAE(
        n_input=n_input,
        n_sample=n_sample,
        n_batch=n_batch,
        n_labels=n_labels,
        n_latent=n_latent,
        n_latent_u=n_latent_u,
        training=False,
    )

    tensors_jax = {
        REGISTRY_KEYS.X_KEY: jnp.array(x_np),
        REGISTRY_KEYS.SAMPLE_KEY: jnp.array(sample_np),
        REGISTRY_KEYS.BATCH_KEY: jnp.array(batch_np),
        REGISTRY_KEYS.LABELS_KEY: jnp.array(label_np),
    }

    key = jax.random.PRNGKey(42)
    keys = jax.random.split(key, 4)
    rngs = {"params": keys[0], "u": keys[1], "dropout": keys[2], "eps": keys[3]}

    variables = jax_mod.init(rngs, tensors_jax)
    jax_params = variables["params"]

    # Define JAX loss function for jax.grad
    def jax_loss_fn(params):
        variables_inner = {"params": params}
        inf_out = jax_mod.apply(
            variables_inner,
            x=jnp.array(x_np),
            sample_index=jnp.array(sample_np),
            use_mean=True,
            method=jax_mod.inference,
            rngs=rngs,
        )
        gen_out = jax_mod.apply(
            variables_inner,
            z=inf_out["z"],
            library=inf_out["library"],
            batch_index=jnp.array(batch_np),
            label_index=jnp.array(label_np),
            method=jax_mod.generative,
            rngs=rngs,
        )
        loss_out = jax_mod.apply(
            variables_inner,
            tensors=tensors_jax,
            inference_outputs=inf_out,
            generative_outputs=gen_out,
            kl_weight=1.0,
            method=jax_mod.loss,
            rngs=rngs,
        )
        return loss_out.loss

    print("\nComputing JAX gradients...")
    jax_grads = jax.grad(jax_loss_fn)(jax_params)

    # ── PyTorch: transfer weights, forward, backward ──
    torch_mod = TorchMRVAE(
        n_input=n_input,
        n_sample=n_sample,
        n_batch=n_batch,
        n_labels=n_labels,
        n_latent=n_latent,
        n_latent_u=n_latent_u,
    )
    torch_mod.eval()  # deterministic (no dropout)
    transfer_all_weights(jax_params, torch_mod)

    tensors_torch = {
        REGISTRY_KEYS.X_KEY: torch.tensor(x_np),
        REGISTRY_KEYS.SAMPLE_KEY: torch.tensor(sample_np),
        REGISTRY_KEYS.BATCH_KEY: torch.tensor(batch_np),
        REGISTRY_KEYS.LABELS_KEY: torch.tensor(label_np),
    }

    print("Computing PyTorch gradients...")
    torch_mod.zero_grad()
    torch_inf = torch_mod.inference(
        x=torch.tensor(x_np),
        sample_index=torch.tensor(sample_np),
        use_mean=True,
    )
    torch_gen = torch_mod.generative(
        z=torch_inf["z"],
        library=torch_inf["library"],
        batch_index=torch.tensor(batch_np),
        label_index=torch.tensor(label_np),
    )
    torch_loss_out = torch_mod.loss(
        tensors=tensors_torch,
        inference_outputs=torch_inf,
        generative_outputs=torch_gen,
        kl_weight=1.0,
    )
    torch_loss_out.loss.backward()

    # ── Compare every gradient ──
    print("\n=== GRADIENT COMPARISON ===")
    compare_all_gradients(jax_grads, torch_mod, ATOL_SIMPLE, ATOL_ATTN)

    print("\n" + "=" * 72)
    print("GRADIENT EQUIVALENCE TEST PASSED")
    print("=" * 72)


if __name__ == "__main__":
    test_full_model_equivalence()
    test_gradient_equivalence()
