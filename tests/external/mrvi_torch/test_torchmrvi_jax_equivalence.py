"""Numerical equivalence tests between JAX and PyTorch MRVI implementations.

Transfers weights from JAX to PyTorch and compares forward pass outputs and
gradients. Requires JAX — skipped automatically if not installed.
"""

from __future__ import annotations

import os

import numpy as np
import pytest
import torch

os.environ.setdefault("JAX_DEFAULT_MATMUL_PRECISION", "float32")
os.environ.setdefault("JAX_PLATFORMS", "cpu")
torch.backends.cuda.matmul.allow_tf32 = False
torch.backends.cudnn.allow_tf32 = False
torch.set_float32_matmul_precision("highest")

from scvi import REGISTRY_KEYS  # noqa: E402
from scvi.external.mrvi_torch._module import TorchMRVAE  # noqa: E402

ATOL_SIMPLE = 1e-5
ATOL_ATTN = 5e-4


# ── weight transfer helpers ──────────────────────────────────────────────────


def _t(arr):
    return torch.tensor(np.array(arr)).float()


def _dense(jax_p, torch_linear):
    torch_linear.weight.data = _t(jax_p["kernel"].T)
    if "bias" in jax_p:
        torch_linear.bias.data = _t(jax_p["bias"])


def _ln(jax_p, torch_ln):
    if "scale" in jax_p:
        torch_ln.weight.data = _t(jax_p["scale"])
    if "bias" in jax_p:
        torch_ln.bias.data = _t(jax_p["bias"])


def _embed(jax_p, torch_emb):
    torch_emb.weight.data = _t(jax_p["embedding"])


def _resnet_block(jp, tb):
    _dense(jp["Dense_0"], tb.fc1)
    _ln(jp["LayerNorm_0"], tb.layer_norm1)
    if tb.fc_match is not None:
        _dense(jp["Dense_1"], tb.fc_match)
        _dense(jp["Dense_2"], tb.fc2)
    else:
        _dense(jp["Dense_1"], tb.fc2)
    _ln(jp["LayerNorm_1"], tb.layer_norm2)


def _mlp(jp, tm):
    for i in range(len(tm.resnet_blocks)):
        _resnet_block(jp[f"ResnetBlock_{i}"], tm.resnet_blocks[i])
    _dense(jp["Dense_0"], tm.fc)


def _normal_dist_nn(jp, tn):
    for i, rb in enumerate(tn.resnet_blocks):
        _resnet_block(jp[f"ResnetBlock_{i}"], rb)
    _dense(jp["Dense_0"], tn.fc_mean)
    _dense(jp["Dense_1"], tn.fc_scale[0])


def _attention_block(jp, tab):
    n_heads, depth, n_channels = tab.n_heads, tab.depth_per_head, tab.n_channels
    tab.query_proj.weight.data = _t(jp["DenseGeneral_0"]["kernel"][:, :, 0].T)
    tab.kv_proj.weight.data = _t(jp["DenseGeneral_1"]["kernel"][:, :, 0].T)
    mha = jp["MultiHeadDotProductAttention_0"]
    for name, proj in [("query", tab.q_proj), ("key", tab.k_proj), ("value", tab.v_proj)]:
        k = np.array(mha[name]["kernel"])
        b = np.array(mha[name]["bias"])
        proj.weight.data = _t(k.reshape(1, n_heads * depth).T)
        proj.bias.data = _t(b.reshape(-1))
    ok = np.array(mha["out"]["kernel"])
    ob = np.array(mha["out"]["bias"])
    tab.out_proj.weight.data = _t(ok.reshape(n_heads * depth, n_channels).T)
    tab.out_proj.bias.data = _t(ob)
    _mlp(jp["MLP_0"], tab.mlp_eps)
    _mlp(jp["MLP_1"], tab.mlp_residual)


def transfer_all_weights(jax_params, torch_module):
    p = jax_params
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

    tqz = torch_module.qz
    _ln(p["qz"]["u_ln"], tqz.layer_norm)
    _embed(p["qz"]["Embed_0"], tqz.embedding)
    _ln(p["qz"]["sample_embed_ln"], tqz.layer_norm_embed)
    _attention_block(p["qz"]["AttentionBlock_0"], tqz.attention_block)
    if tqz.fc is not None:
        _dense(p["qz"]["Dense_0"], tqz.fc)

    tpx = torch_module.px
    _ln(p["px"]["u_ln"], tpx.layer_norm)
    _embed(p["px"]["Embed_0"], tpx.batch_embedding)
    _ln(p["px"]["batch_embed_ln"], tpx.layer_norm_batch_embed)
    _attention_block(p["px"]["AttentionBlock_0"], tpx.attention_block)
    _dense(p["px"]["Dense_0"], tpx.fc)
    tpx.px_r.data = _t(p["px"]["px_r"])

    if torch_module.u_prior_mixture:
        torch_module.u_prior_logits.data = _t(p["u_prior_logits"])
        torch_module.u_prior_means.data = _t(p["u_prior_means"])
        torch_module.u_prior_scales.data = _t(p["u_prior_scales"])

    if hasattr(torch_module, "pz_scale"):
        if isinstance(torch_module.pz_scale, torch.nn.Parameter):
            torch_module.pz_scale.data = _t(p["pz_scale"])
        else:
            pz_val = (
                np.array(p["pz_scale"]) if "pz_scale" in p else np.zeros(torch_module.n_latent)
            )
            torch_module.pz_scale.copy_(_t(pz_val))


# ── comparison helpers ───────────────────────────────────────────────────────


def _compare(name, jax_val, torch_val, atol):
    j = np.array(jax_val)
    t = (
        torch_val.detach().cpu().numpy()
        if isinstance(torch_val, torch.Tensor)
        else np.array(torch_val)
    )
    diff = np.max(np.abs(j - t))
    assert diff < atol, f"{name}: max_diff={diff:.2e} exceeds atol={atol:.0e}"


def _compare_grad(name, jax_grad_val, torch_param, atol):
    j = np.array(jax_grad_val)
    assert torch_param.grad is not None, f"{name}: PyTorch grad is None"
    t = torch_param.grad.detach().cpu().numpy()
    diff = np.max(np.abs(j - t))
    assert diff < atol, f"grad {name}: max_diff={diff:.2e} exceeds atol={atol:.0e}"


def _compare_dense_grad(name, jg, tl, atol):
    _compare_grad(f"{name}.w", np.array(jg["kernel"]).T, tl.weight, atol)
    if "bias" in jg and tl.bias is not None:
        _compare_grad(f"{name}.b", jg["bias"], tl.bias, atol)


def _compare_ln_grad(name, jg, tl, atol):
    if "scale" in jg:
        _compare_grad(f"{name}.w", jg["scale"], tl.weight, atol)
    if "bias" in jg:
        _compare_grad(f"{name}.b", jg["bias"], tl.bias, atol)


def _compare_embed_grad(name, jg, te, atol):
    _compare_grad(f"{name}.w", jg["embedding"], te.weight, atol)


def _compare_resnet_block_grad(name, jg, tb, atol):
    _compare_dense_grad(f"{name}.fc1", jg["Dense_0"], tb.fc1, atol)
    _compare_ln_grad(f"{name}.ln1", jg["LayerNorm_0"], tb.layer_norm1, atol)
    if tb.fc_match is not None:
        _compare_dense_grad(f"{name}.fc_match", jg["Dense_1"], tb.fc_match, atol)
        _compare_dense_grad(f"{name}.fc2", jg["Dense_2"], tb.fc2, atol)
    else:
        _compare_dense_grad(f"{name}.fc2", jg["Dense_1"], tb.fc2, atol)
    _compare_ln_grad(f"{name}.ln2", jg["LayerNorm_1"], tb.layer_norm2, atol)


def _compare_mlp_grad(name, jg, tm, atol):
    for i in range(len(tm.resnet_blocks)):
        _compare_resnet_block_grad(
            f"{name}.rb{i}", jg[f"ResnetBlock_{i}"], tm.resnet_blocks[i], atol
        )
    _compare_dense_grad(f"{name}.fc", jg["Dense_0"], tm.fc, atol)


def _compare_attention_block_grad(name, jg, tab, atol):
    n_heads, depth, n_channels = tab.n_heads, tab.depth_per_head, tab.n_channels
    _compare_grad(
        f"{name}.qproj.w",
        np.array(jg["DenseGeneral_0"]["kernel"])[:, :, 0].T,
        tab.query_proj.weight,
        atol,
    )
    _compare_grad(
        f"{name}.kvproj.w",
        np.array(jg["DenseGeneral_1"]["kernel"])[:, :, 0].T,
        tab.kv_proj.weight,
        atol,
    )
    mha_g = jg["MultiHeadDotProductAttention_0"]
    for qkv, proj in [("query", tab.q_proj), ("key", tab.k_proj), ("value", tab.v_proj)]:
        _compare_grad(
            f"{name}.{qkv}.w",
            np.array(mha_g[qkv]["kernel"]).reshape(1, n_heads * depth).T,
            proj.weight,
            atol,
        )
        _compare_grad(f"{name}.{qkv}.b", np.array(mha_g[qkv]["bias"]).reshape(-1), proj.bias, atol)
    _compare_grad(
        f"{name}.out.w",
        np.array(mha_g["out"]["kernel"]).reshape(n_heads * depth, n_channels).T,
        tab.out_proj.weight,
        atol,
    )
    _compare_grad(f"{name}.out.b", np.array(mha_g["out"]["bias"]), tab.out_proj.bias, atol)
    _compare_mlp_grad(f"{name}.mlp_eps", jg["MLP_0"], tab.mlp_eps, atol)
    _compare_mlp_grad(f"{name}.mlp_res", jg["MLP_1"], tab.mlp_residual, atol)


def compare_all_gradients(jax_grads, torch_mod, atol_s, atol_a):
    g = jax_grads
    tqu = torch_mod.qu
    _compare_dense_grad("qu.fc1", g["qu"]["Dense_0"], tqu.fc1, atol_s)
    _compare_dense_grad("qu.fc2", g["qu"]["Dense_1"], tqu.fc2, atol_s)
    for i, cn in enumerate([tqu.conditional_norm1, tqu.conditional_norm2]):
        _compare_embed_grad(
            f"qu.cn{i}.gamma",
            g["qu"][f"ConditionalNormalization_{i}"]["gamma_conditional"],
            cn.gamma_embedding,
            atol_s,
        )
        _compare_embed_grad(
            f"qu.cn{i}.beta",
            g["qu"][f"ConditionalNormalization_{i}"]["beta_conditional"],
            cn.beta_embedding,
            atol_s,
        )
    _compare_embed_grad("qu.sample_embed", g["qu"]["Embed_0"], tqu.sample_embed, atol_s)
    nn_g = g["qu"]["NormalDistOutputNN_0"]
    for i, rb in enumerate(tqu.output_nn.resnet_blocks):
        _compare_resnet_block_grad(f"qu.nn.rb{i}", nn_g[f"ResnetBlock_{i}"], rb, atol_s)
    _compare_dense_grad("qu.nn.fc_mean", nn_g["Dense_0"], tqu.output_nn.fc_mean, atol_s)
    _compare_dense_grad("qu.nn.fc_scale", nn_g["Dense_1"], tqu.output_nn.fc_scale[0], atol_s)

    tqz = torch_mod.qz
    _compare_ln_grad("qz.u_ln", g["qz"]["u_ln"], tqz.layer_norm, atol_a)
    _compare_embed_grad("qz.embed", g["qz"]["Embed_0"], tqz.embedding, atol_a)
    _compare_ln_grad("qz.embed_ln", g["qz"]["sample_embed_ln"], tqz.layer_norm_embed, atol_a)
    _compare_attention_block_grad(
        "qz.attn", g["qz"]["AttentionBlock_0"], tqz.attention_block, atol_a
    )
    if tqz.fc is not None:
        _compare_dense_grad("qz.fc", g["qz"]["Dense_0"], tqz.fc, atol_a)

    tpx = torch_mod.px
    _compare_ln_grad("px.u_ln", g["px"]["u_ln"], tpx.layer_norm, atol_a)
    _compare_embed_grad("px.batch_embed", g["px"]["Embed_0"], tpx.batch_embedding, atol_a)
    _compare_ln_grad("px.batch_ln", g["px"]["batch_embed_ln"], tpx.layer_norm_batch_embed, atol_a)
    _compare_attention_block_grad(
        "px.attn", g["px"]["AttentionBlock_0"], tpx.attention_block, atol_a
    )
    _compare_dense_grad("px.fc", g["px"]["Dense_0"], tpx.fc, atol_a)
    _compare_grad("px.px_r", g["px"]["px_r"], tpx.px_r, atol_a)

    if torch_mod.u_prior_mixture:
        _compare_grad("u_prior_logits", g["u_prior_logits"], torch_mod.u_prior_logits, atol_s)
        _compare_grad("u_prior_means", g["u_prior_means"], torch_mod.u_prior_means, atol_s)
        _compare_grad("u_prior_scales", g["u_prior_scales"], torch_mod.u_prior_scales, atol_s)


# ── shared test data ─────────────────────────────────────────────────────────

_N_INPUT, _N_SAMPLE, _N_BATCH, _N_LABELS = 100, 5, 2, 1
_N_LATENT, _N_LATENT_U = 30, 10
_BS = 16


def _make_test_data():
    np.random.seed(0)
    return {
        "x": (np.abs(np.random.randn(_BS, _N_INPUT)) + 0.1).astype(np.float32),
        "sample": np.random.randint(0, _N_SAMPLE, (_BS, 1)).astype(np.float32),
        "batch": np.random.randint(0, _N_BATCH, (_BS, 1)).astype(np.float32),
        "label": np.zeros((_BS, 1), dtype=np.float32),
    }


def _init_jax_module(data):
    jax = pytest.importorskip("jax")
    jnp = pytest.importorskip("jax.numpy")
    JaxMRVAE = pytest.importorskip("scvi.external.mrvi_jax._module").JaxMRVAE

    jax.config.update("jax_default_matmul_precision", "float32")

    jax_mod = JaxMRVAE(
        n_input=_N_INPUT,
        n_sample=_N_SAMPLE,
        n_batch=_N_BATCH,
        n_labels=_N_LABELS,
        n_latent=_N_LATENT,
        n_latent_u=_N_LATENT_U,
        training=False,
    )
    tensors = {
        REGISTRY_KEYS.X_KEY: jnp.array(data["x"]),
        REGISTRY_KEYS.SAMPLE_KEY: jnp.array(data["sample"]),
        REGISTRY_KEYS.BATCH_KEY: jnp.array(data["batch"]),
        REGISTRY_KEYS.LABELS_KEY: jnp.array(data["label"]),
    }
    key = jax.random.PRNGKey(42)
    keys = jax.random.split(key, 4)
    rngs = {"params": keys[0], "u": keys[1], "dropout": keys[2], "eps": keys[3]}
    variables = jax_mod.init(rngs, tensors)
    return jax_mod, variables, tensors, rngs


def _init_torch_module(jax_params):
    torch_mod = TorchMRVAE(
        n_input=_N_INPUT,
        n_sample=_N_SAMPLE,
        n_batch=_N_BATCH,
        n_labels=_N_LABELS,
        n_latent=_N_LATENT,
        n_latent_u=_N_LATENT_U,
    )
    torch_mod.eval()
    transfer_all_weights(jax_params, torch_mod)
    return torch_mod


# ── tests ────────────────────────────────────────────────────────────────────


@pytest.mark.jax
def test_forward_pass_equivalence():
    """Init JAX, transfer weights to PyTorch, compare all inference/generative/loss outputs."""
    pytest.importorskip("jax")
    jnp = pytest.importorskip("jax.numpy")

    data = _make_test_data()
    jax_mod, variables, tensors_jax, rngs = _init_jax_module(data)
    torch_mod = _init_torch_module(variables["params"])

    # Inference
    jax_inf = jax_mod.apply(
        variables,
        x=jnp.array(data["x"]),
        sample_index=jnp.array(data["sample"]),
        use_mean=True,
        method=jax_mod.inference,
        rngs=rngs,
    )
    with torch.no_grad():
        torch_inf = torch_mod.inference(
            x=torch.tensor(data["x"]), sample_index=torch.tensor(data["sample"]), use_mean=True
        )

    for key in ("u", "z", "z_base", "library"):
        atol = ATOL_SIMPLE if key in ("u", "library") else ATOL_ATTN
        _compare(key, jax_inf[key], torch_inf[key], atol)
    _compare("qu.mean", jax_inf["qu"].mean, torch_inf["qu"].mean, ATOL_SIMPLE)
    _compare("qu.scale", jax_inf["qu"].scale, torch_inf["qu"].scale, ATOL_SIMPLE)

    # Generative
    jax_gen = jax_mod.apply(
        variables,
        z=jax_inf["z"],
        library=jax_inf["library"],
        batch_index=jnp.array(data["batch"]),
        label_index=jnp.array(data["label"]),
        method=jax_mod.generative,
        rngs=rngs,
    )
    with torch.no_grad():
        torch_gen = torch_mod.generative(
            z=torch_inf["z"],
            library=torch_inf["library"],
            batch_index=torch.tensor(data["batch"]),
            label_index=torch.tensor(data["label"]),
        )
    _compare("h", jax_gen["h"], torch_gen["h"], ATOL_ATTN)
    _compare("px.mean", jax_gen["px"].mean, torch_gen["px"].mean, ATOL_ATTN)

    # Loss
    tensors_torch = {
        REGISTRY_KEYS.X_KEY: torch.tensor(data["x"]),
        REGISTRY_KEYS.SAMPLE_KEY: torch.tensor(data["sample"]),
        REGISTRY_KEYS.BATCH_KEY: torch.tensor(data["batch"]),
        REGISTRY_KEYS.LABELS_KEY: torch.tensor(data["label"]),
    }

    jax_loss = jax_mod.apply(
        variables,
        tensors=tensors_jax,
        inference_outputs=jax_inf,
        generative_outputs=jax_gen,
        kl_weight=1.0,
        method=jax_mod.loss,
        rngs=rngs,
    )
    with torch.no_grad():
        torch_loss = torch_mod.loss(
            tensors=tensors_torch,
            inference_outputs=torch_inf,
            generative_outputs=torch_gen,
            kl_weight=1.0,
        )

    def _extract(val):
        return sum(val.values()) if isinstance(val, dict) else val

    _compare(
        "recon_loss",
        _extract(jax_loss.reconstruction_loss),
        _extract(torch_loss.reconstruction_loss),
        ATOL_ATTN,
    )
    _compare("kl_local", _extract(jax_loss.kl_local), _extract(torch_loss.kl_local), ATOL_ATTN)
    _compare("total_loss", jax_loss.loss, torch_loss.loss, ATOL_ATTN)


@pytest.mark.jax
def test_gradient_equivalence():
    """Compare jax.grad vs loss.backward() for all parameters."""
    jax = pytest.importorskip("jax")
    jnp = pytest.importorskip("jax.numpy")

    data = _make_test_data()
    jax_mod, variables, tensors_jax, rngs = _init_jax_module(data)

    def jax_loss_fn(params):
        v = {"params": params}
        inf = jax_mod.apply(
            v,
            x=jnp.array(data["x"]),
            sample_index=jnp.array(data["sample"]),
            use_mean=True,
            method=jax_mod.inference,
            rngs=rngs,
        )
        gen = jax_mod.apply(
            v,
            z=inf["z"],
            library=inf["library"],
            batch_index=jnp.array(data["batch"]),
            label_index=jnp.array(data["label"]),
            method=jax_mod.generative,
            rngs=rngs,
        )
        loss = jax_mod.apply(
            v,
            tensors=tensors_jax,
            inference_outputs=inf,
            generative_outputs=gen,
            kl_weight=1.0,
            method=jax_mod.loss,
            rngs=rngs,
        )
        return loss.loss

    jax_grads = jax.grad(jax_loss_fn)(variables["params"])

    torch_mod = _init_torch_module(variables["params"])
    torch_mod.zero_grad()
    torch_inf = torch_mod.inference(
        x=torch.tensor(data["x"]), sample_index=torch.tensor(data["sample"]), use_mean=True
    )
    torch_gen = torch_mod.generative(
        z=torch_inf["z"],
        library=torch_inf["library"],
        batch_index=torch.tensor(data["batch"]),
        label_index=torch.tensor(data["label"]),
    )
    tensors_torch = {
        REGISTRY_KEYS.X_KEY: torch.tensor(data["x"]),
        REGISTRY_KEYS.SAMPLE_KEY: torch.tensor(data["sample"]),
        REGISTRY_KEYS.BATCH_KEY: torch.tensor(data["batch"]),
        REGISTRY_KEYS.LABELS_KEY: torch.tensor(data["label"]),
    }
    torch_mod.loss(
        tensors=tensors_torch,
        inference_outputs=torch_inf,
        generative_outputs=torch_gen,
        kl_weight=1.0,
    ).loss.backward()

    compare_all_gradients(jax_grads, torch_mod, ATOL_SIMPLE, ATOL_ATTN)
