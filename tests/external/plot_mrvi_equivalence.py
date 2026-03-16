"""Visualize numerical equivalence between JAX and PyTorch MRVI implementations.

Produces two figures:
  1. Parameter update equivalence across multiple gradient steps
  2. Differential expression result equivalence on synthetic data

Both implementations start from identical weights, use the same data,
and results are compared element-wise.
"""

from __future__ import annotations

import os

os.environ["JAX_DEFAULT_MATMUL_PRECISION"] = "float32"

import warnings

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import torch

torch.backends.cuda.matmul.allow_tf32 = False
torch.backends.cudnn.allow_tf32 = False
torch.set_float32_matmul_precision("highest")

import jax  # noqa: E402
import jax.numpy as jnp  # noqa: E402
import optax  # noqa: E402

jax.config.update("jax_default_matmul_precision", "float32")

# Reuse weight transfer from the equivalence test
import sys  # noqa: E402

from scvi import REGISTRY_KEYS  # noqa: E402
from scvi.external.mrvi_jax._module import JaxMRVAE  # noqa: E402
from scvi.external.mrvi_torch._module import TorchMRVAE  # noqa: E402

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from test_mrvi_full_model_equivalence import (  # noqa: E402
    transfer_all_weights,
)

warnings.filterwarnings("ignore")

# ── Helpers ──────────────────────────────────────────────────────────────────


def _flat_params_jax(params):
    """Flatten JAX param tree to a dict of path→numpy array."""
    out = {}
    for path, val in jax.tree_util.tree_leaves_with_path(params):
        key = "/".join(str(getattr(k, "key", k)) for k in path)
        out[key] = np.array(val)
    return out


def _flat_params_torch(module):
    """Flatten PyTorch named parameters to a dict of name→numpy array."""
    return {name: p.detach().cpu().numpy().copy() for name, p in module.named_parameters()}


def _collect_torch_grads(module):
    """Collect PyTorch gradients as a dict."""
    return {
        name: p.grad.detach().cpu().numpy().copy()
        for name, p in module.named_parameters()
        if p.grad is not None
    }


# ── 1. Parameter update equivalence ─────────────────────────────────────────


def run_parameter_update_comparison(n_steps=10):
    """Run gradient steps on both backends, track parameter diffs at each step."""
    n_input, n_sample, n_batch, n_labels = 50, 3, 2, 1
    n_latent, n_latent_u = 10, 5
    bs = 32
    lr = 1e-3

    np.random.seed(0)
    x_np = (np.abs(np.random.randn(bs, n_input)) + 0.1).astype(np.float32)
    sample_np = np.random.randint(0, n_sample, (bs, 1)).astype(np.float32)
    batch_np = np.random.randint(0, n_batch, (bs, 1)).astype(np.float32)
    label_np = np.zeros((bs, 1), dtype=np.float32)

    tensors_jax = {
        REGISTRY_KEYS.X_KEY: jnp.array(x_np),
        REGISTRY_KEYS.SAMPLE_KEY: jnp.array(sample_np),
        REGISTRY_KEYS.BATCH_KEY: jnp.array(batch_np),
        REGISTRY_KEYS.LABELS_KEY: jnp.array(label_np),
    }

    # ── JAX setup ──
    jax_mod = JaxMRVAE(
        n_input=n_input,
        n_sample=n_sample,
        n_batch=n_batch,
        n_labels=n_labels,
        n_latent=n_latent,
        n_latent_u=n_latent_u,
        training=False,
    )
    key = jax.random.PRNGKey(42)
    keys = jax.random.split(key, 4)
    rngs = {"params": keys[0], "u": keys[1], "dropout": keys[2], "eps": keys[3]}
    variables = jax_mod.init(rngs, tensors_jax)
    jax_params = variables["params"]

    # JAX optimizer (SGD to keep it simple and deterministic)
    jax_opt = optax.sgd(lr)
    jax_opt_state = jax_opt.init(jax_params)

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

    # ── PyTorch setup ──
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
    torch_opt = torch.optim.SGD(torch_mod.parameters(), lr=lr)

    tensors_torch = {
        REGISTRY_KEYS.X_KEY: torch.tensor(x_np),
        REGISTRY_KEYS.SAMPLE_KEY: torch.tensor(sample_np),
        REGISTRY_KEYS.BATCH_KEY: torch.tensor(batch_np),
        REGISTRY_KEYS.LABELS_KEY: torch.tensor(label_np),
    }

    # ── Training loop ──
    jax_losses = []
    torch_losses = []
    param_max_diffs = []  # max diff across ALL params at each step
    param_mean_diffs = []
    grad_max_diffs = []

    for step in range(n_steps):
        # ── JAX step ──
        jax_loss, jax_grads = jax.value_and_grad(jax_loss_fn)(jax_params)
        updates, jax_opt_state = jax_opt.update(jax_grads, jax_opt_state)
        jax_params = optax.apply_updates(jax_params, updates)
        jax_losses.append(float(jax_loss))

        # ── PyTorch step ──
        torch_opt.zero_grad()
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
        torch_opt.step()
        torch_losses.append(float(torch_loss_out.loss.item()))

        # ── Compare gradients (before param update applied in torch) ──
        # For JAX, grads were computed before update. For Torch, .grad is populated before step().
        # Since we already called step(), we compare the grads from this iteration.
        _flat_params_jax(jax_grads)
        torch_grad_flat = _collect_torch_grads(torch_mod)
        # Grad comparison is approximate because param names differ between JAX/PyTorch
        # Use aggregate statistics instead
        all_grad_diffs = []
        for _tname, _tgrad in torch_grad_flat.items():
            all_grad_diffs.append(0.0)  # placeholder — we compare aggregate below
        # Simpler: compare loss values directly
        grad_max_diffs.append(abs(jax_losses[-1] - torch_losses[-1]))

        # ── Compare parameters after update ──
        # Transfer updated JAX params to a fresh torch module for comparison
        _flat_params_jax(jax_params)
        torch_flat = _flat_params_torch(torch_mod)
        all_diffs = []
        for _tname, _tval in torch_flat.items():
            all_diffs.append(0.0)  # individual diffs computed below

        # Compute aggregate param diff by re-transferring JAX params to a reference
        # and comparing against the trained torch module
        ref_mod = TorchMRVAE(
            n_input=n_input,
            n_sample=n_sample,
            n_batch=n_batch,
            n_labels=n_labels,
            n_latent=n_latent,
            n_latent_u=n_latent_u,
        )
        ref_mod.eval()
        transfer_all_weights(jax_params, ref_mod)
        ref_flat = _flat_params_torch(ref_mod)
        diffs = []
        for name in torch_flat:
            if name in ref_flat and torch_flat[name].shape == ref_flat[name].shape:
                diffs.append(np.max(np.abs(torch_flat[name] - ref_flat[name])))
        param_max_diffs.append(max(diffs) if diffs else 0.0)
        param_mean_diffs.append(np.mean(diffs) if diffs else 0.0)

        print(
            f"Step {step:3d}  |  JAX loss={jax_losses[-1]:.6f}  "
            f"Torch loss={torch_losses[-1]:.6f}  "
            f"|loss diff|={abs(jax_losses[-1] - torch_losses[-1]):.2e}  "
            f"param max_diff={param_max_diffs[-1]:.2e}"
        )

    return {
        "jax_losses": jax_losses,
        "torch_losses": torch_losses,
        "param_max_diffs": param_max_diffs,
        "param_mean_diffs": param_mean_diffs,
        "n_steps": n_steps,
    }


# ── 2. DE equivalence ───────────────────────────────────────────────────────


def run_de_comparison():
    """Train JAX model, transfer weights to PyTorch, run DE on both, compare.

    Both models share identical weights. DE involves stochastic MC sampling with
    different RNG systems, so beta/effect_size/lfc won't match exactly, but should
    be highly correlated since they use the same model and same data.
    """
    import scvi
    from scvi.external import MRVI

    np.random.seed(42)
    adata = scvi.data.synthetic_iid(n_genes=200)
    n_obs = adata.n_obs
    adata.obs["sample"] = np.random.choice(8, n_obs).astype(str)
    meta_per_sample = {str(i): int(i >= 4) for i in range(8)}
    adata.obs["meta1"] = adata.obs["sample"].map(meta_per_sample).astype(int)
    adata.obs["meta1_cat"] = "CAT_" + adata.obs["meta1"].astype(str)
    adata.obs["meta1_cat"] = adata.obs["meta1_cat"].astype("category")

    # Inject a real sample-level covariate effect: upregulate genes 0-19 for CAT_1
    import scipy.sparse as sp

    x = adata.X.toarray() if sp.issparse(adata.X) else adata.X.copy()
    for sample_id in range(8):
        mask = adata.obs["sample"] == str(sample_id)
        if meta_per_sample[str(sample_id)] == 1:
            x[mask, :20] = x[mask, :20] * 3.0 + 5.0
    adata.X = x.astype(np.float32)

    # ── Train JAX model ──
    print("\nTraining JAX model...")
    MRVI.setup_anndata(adata, sample_key="sample", batch_key="batch", backend="jax")
    jax_model = MRVI(adata)
    jax_model.train(max_epochs=20, batch_size=256, train_size=0.9)
    jax_model.update_sample_info(adata)

    # ── Create PyTorch model with identical weights ──
    print("Transferring JAX weights to PyTorch model...")
    MRVI.setup_anndata(adata, sample_key="sample", batch_key="batch", backend="torch")
    torch_model = MRVI(adata)
    jax_params = jax_model.module.params
    transfer_all_weights(jax_params, torch_model.module)
    torch_model.module.eval()
    torch_model.is_trained_ = True
    torch_model.update_sample_info(adata)

    # ── Compute local sample distances on both (deterministic, use_mean=True) ──
    dist_kwargs = {"batch_size": 400, "use_mean": True, "use_vmap": False}
    print("Computing JAX local sample distances...")
    jax_dists = jax_model.get_local_sample_distances(**dist_kwargs)
    print("Computing PyTorch local sample distances...")
    torch_dists = torch_model.get_local_sample_distances(**dist_kwargs)

    # ── Run DE on both ──
    # DE uses use_mean=True for the counterfactual representations (deterministic),
    # but internally uses MC sampling for LFC computation. We use high mc_samples
    # to minimize variance.
    de_kwargs = {
        "sample_cov_keys": ["meta1_cat"],
        "store_lfc": True,
        "mc_samples": 500,
        "batch_size": 400,
        "use_vmap": False,
    }
    print("Running JAX DE...")
    jax_de = jax_model.differential_expression(**de_kwargs)

    print("Running PyTorch DE...")
    torch_de = torch_model.differential_expression(**de_kwargs)

    # ── Also compute deterministic local sample representations ──
    # These are the inputs to the DE beta computation and are fully deterministic.
    print("Computing JAX local sample representations (use_mean=True)...")
    jax_reps = jax_model.get_local_sample_representation(batch_size=400, use_mean=True)
    print("Computing PyTorch local sample representations (use_mean=True)...")
    torch_reps = torch_model.get_local_sample_representation(batch_size=400, use_mean=True)

    return jax_de, torch_de, jax_dists, torch_dists, jax_reps, torch_reps


# ── Plotting ─────────────────────────────────────────────────────────────────


def plot_parameter_updates(results, save_path):
    """Plot parameter update equivalence with before/after comparison."""
    # Load old (pre-fix) results if available
    old_results_path = "/tmp/mrvi_old_results.npz"
    has_old = os.path.exists(old_results_path)
    if has_old:
        old = np.load(old_results_path)
        old_loss_diffs = np.abs(old["jax_losses"] - old["torch_losses"])

    fig, axes = plt.subplots(1, 3, figsize=(16, 4.5))
    steps = np.arange(results["n_steps"])

    # Panel 1: Loss curves (after fix)
    ax = axes[0]
    ax.plot(steps, results["jax_losses"], "o-", label="JAX", markersize=4, color="#2196F3")
    ax.plot(steps, results["torch_losses"], "s--", label="PyTorch", markersize=4, color="#FF5722")
    ax.set_xlabel("Gradient step")
    ax.set_ylabel("Loss")
    ax.set_title("Training loss (after fix)")
    ax.legend()
    ax.grid(True, alpha=0.3)

    # Panel 2: Loss difference — before vs after
    ax = axes[1]
    loss_diffs = np.abs(np.array(results["jax_losses"]) - np.array(results["torch_losses"]))
    if has_old:
        ax.semilogy(
            steps,
            old_loss_diffs[: len(steps)],
            "x--",
            color="#BDBDBD",
            markersize=6,
            linewidth=2,
            label="Before fix",
            zorder=1,
        )
    ax.semilogy(
        steps,
        loss_diffs,
        "o-",
        color="#4CAF50",
        markersize=4,
        label="After fix",
        zorder=2,
    )
    ax.set_xlabel("Gradient step")
    ax.set_ylabel("|JAX loss - PyTorch loss|")
    ax.set_title("Loss difference (log scale)")
    ax.grid(True, alpha=0.3)
    ax.legend()

    # Panel 3: Parameter max diff (after fix only — old code can't transfer weights)
    ax = axes[2]
    ax.semilogy(
        steps, results["param_max_diffs"], "o-", label="Max diff", color="#9C27B0", markersize=4
    )
    ax.semilogy(
        steps, results["param_mean_diffs"], "s--", label="Mean diff", color="#FF9800", markersize=4
    )
    ax.set_xlabel("Gradient step")
    ax.set_ylabel("Parameter difference")
    ax.set_title("Parameter divergence (log scale)")
    ax.legend()
    ax.grid(True, alpha=0.3)
    ax.axhline(y=1e-4, color="gray", linestyle=":", alpha=0.5)

    fig.suptitle("JAX vs PyTorch MRVI: Parameter Update Equivalence", fontsize=13, y=1.02)
    fig.tight_layout()
    fig.savefig(save_path, dpi=150, bbox_inches="tight")
    print(f"Saved parameter update plot to {save_path}")
    plt.close(fig)


def plot_de_comparison(jax_de, torch_de, jax_reps, torch_reps, save_path):
    """Plot DE comparison: JAX vs PyTorch, same weights.

    Includes a deterministic representation panel to prove the forward pass
    matches, alongside stochastic DE results (which have MC sampling noise).
    """
    fig, axes = plt.subplots(2, 3, figsize=(16, 9))

    # ── Row 1, Panel 1: Deterministic local sample representations ──
    ax = axes[0, 0]
    j = jax_reps.values.flatten()
    t = torch_reps.values.flatten()
    if len(j) > 50000:
        idx = np.random.choice(len(j), 50000, replace=False)
        jp, tp = j[idx], t[idx]
    else:
        jp, tp = j, t
    ax.scatter(jp, tp, s=2, alpha=0.2, color="#2196F3", edgecolors="none")
    lim = max(np.abs(j).max(), np.abs(t).max()) * 1.1
    ax.plot([-lim, lim], [-lim, lim], "k--", lw=0.8, alpha=0.5)
    ax.set(xlabel="JAX", ylabel="PyTorch", xlim=(-lim, lim), ylim=(-lim, lim))
    ax.set_aspect("equal")
    corr_rep = np.corrcoef(j, t)[0, 1]
    max_diff_rep = np.max(np.abs(j - t))
    ax.set_title(
        f"Local sample representations\n(deterministic, use_mean=True)\n"
        f"r={corr_rep:.6f}, max_diff={max_diff_rep:.2e}"
    )
    ax.grid(True, alpha=0.3)

    # ── Row 1, Panel 2: Beta coefficients ──
    ax = axes[0, 1]
    j, t = jax_de["beta"].values.flatten(), torch_de["beta"].values.flatten()
    ax.scatter(j, t, s=8, alpha=0.5, color="#4CAF50", edgecolors="none")
    lim = max(np.abs(j).max(), np.abs(t).max()) * 1.1
    ax.plot([-lim, lim], [-lim, lim], "k--", lw=0.8, alpha=0.5)
    ax.set(xlabel="JAX beta", ylabel="PyTorch beta", xlim=(-lim, lim), ylim=(-lim, lim))
    ax.set_aspect("equal")
    ax.set_title(f"DE beta coefficients\nr={np.corrcoef(j, t)[0, 1]:.6f}")
    ax.grid(True, alpha=0.3)

    # ── Row 1, Panel 3: Effect sizes ──
    ax = axes[0, 2]
    j, t = jax_de["effect_size"].values.flatten(), torch_de["effect_size"].values.flatten()
    ax.scatter(j, t, s=8, alpha=0.5, color="#FF9800", edgecolors="none")
    lim = max(j.max(), t.max()) * 1.1
    ax.plot([0, lim], [0, lim], "k--", lw=0.8, alpha=0.5)
    ax.set(xlabel="JAX effect size", ylabel="PyTorch effect size")
    ax.set_title(f"DE effect sizes\nr={np.corrcoef(j, t)[0, 1]:.6f}")
    ax.grid(True, alpha=0.3)

    # ── Row 2, Panel 1: LFC ──
    ax = axes[1, 0]
    if "lfc" in jax_de and "lfc" in torch_de:
        j, t = jax_de["lfc"].values.flatten(), torch_de["lfc"].values.flatten()
        m = np.isfinite(j) & np.isfinite(t)
        j, t = j[m], t[m]
        ax.scatter(j, t, s=4, alpha=0.3, color="#FF5722", edgecolors="none")
        lim = max(np.abs(j).max(), np.abs(t).max()) * 1.1
        ax.plot([-lim, lim], [-lim, lim], "k--", lw=0.8, alpha=0.5)
        ax.set(xlabel="JAX LFC", ylabel="PyTorch LFC", xlim=(-lim, lim), ylim=(-lim, lim))
        ax.set_aspect("equal")
        ax.set_title(f"DE log fold changes\nr={np.corrcoef(j, t)[0, 1]:.6f}")
    ax.grid(True, alpha=0.3)

    # ── Row 2, Panel 2: P-values ──
    ax = axes[1, 1]
    j, t = jax_de["pvalue"].values.flatten(), torch_de["pvalue"].values.flatten()
    ax.scatter(j, t, s=8, alpha=0.5, color="#9C27B0", edgecolors="none")
    ax.plot([0, 1], [0, 1], "k--", lw=0.8, alpha=0.5)
    ax.set(xlabel="JAX p-value", ylabel="PyTorch p-value", xlim=(0, 1), ylim=(0, 1))
    ax.set_aspect("equal")
    pv_corr = np.corrcoef(j, t)[0, 1]
    pv_max_diff = np.max(np.abs(j - t))
    pv_label = f"r={pv_corr:.6f}" if np.isfinite(pv_corr) else f"max_diff={pv_max_diff:.2e}"
    ax.set_title(f"DE p-values\n{pv_label}")
    ax.grid(True, alpha=0.3)

    # ── Row 2, Panel 3: padj ──
    ax = axes[1, 2]
    j, t = jax_de["padj"].values.flatten(), torch_de["padj"].values.flatten()
    ax.scatter(j, t, s=8, alpha=0.5, color="#009688", edgecolors="none")
    ax.plot([0, 1], [0, 1], "k--", lw=0.8, alpha=0.5)
    ax.set(xlabel="JAX padj", ylabel="PyTorch padj", xlim=(0, 1), ylim=(0, 1))
    ax.set_aspect("equal")
    padj_corr = np.corrcoef(j, t)[0, 1]
    padj_max_diff = np.max(np.abs(j - t))
    padj_label = (
        f"r={padj_corr:.6f}" if np.isfinite(padj_corr) else f"max_diff={padj_max_diff:.2e}"
    )
    ax.set_title(f"DE adjusted p-values\n{padj_label}")
    ax.grid(True, alpha=0.3)

    fig.suptitle(
        "JAX vs PyTorch MRVI: Differential Expression (same weights, mc_samples=500)",
        fontsize=13,
        y=1.02,
    )
    fig.tight_layout()
    fig.savefig(save_path, dpi=150, bbox_inches="tight")
    print(f"Saved DE comparison plot to {save_path}")
    plt.close(fig)


def plot_distance_comparison(jax_dists, torch_dists, save_path):
    """Plot local sample distance matrix comparison between JAX and PyTorch."""
    jax_d = jax_dists["cell"].values  # (n_cells, n_samples, n_samples)
    torch_d = torch_dists["cell"].values

    fig, axes = plt.subplots(1, 3, figsize=(16, 4.5))

    # Panel 1: Scatter of all distance matrix elements
    ax = axes[0]
    jax_flat = jax_d.flatten()
    torch_flat = torch_d.flatten()
    # Subsample for plotting if too many points
    n_pts = len(jax_flat)
    if n_pts > 50000:
        idx = np.random.choice(n_pts, 50000, replace=False)
        jax_flat_plot, torch_flat_plot = jax_flat[idx], torch_flat[idx]
    else:
        jax_flat_plot, torch_flat_plot = jax_flat, torch_flat
    ax.scatter(jax_flat_plot, torch_flat_plot, s=2, alpha=0.2, color="#2196F3", edgecolors="none")
    lim = max(jax_flat.max(), torch_flat.max()) * 1.05
    ax.plot([0, lim], [0, lim], "k--", linewidth=0.8, alpha=0.5)
    corr = np.corrcoef(jax_flat, torch_flat)[0, 1]
    max_diff = np.max(np.abs(jax_flat - torch_flat))
    ax.set_xlabel("JAX distances")
    ax.set_ylabel("PyTorch distances")
    ax.set_title(f"All pairwise distances\n(corr={corr:.6f}, max_diff={max_diff:.2e})")
    ax.set_xlim(0, lim)
    ax.set_ylim(0, lim)
    ax.set_aspect("equal")
    ax.grid(True, alpha=0.3)

    # Panel 2: Mean distance matrix (averaged over cells)
    ax = axes[1]
    jax_mean = jax_d.mean(axis=0)  # (n_samples, n_samples)
    torch_mean = torch_d.mean(axis=0)
    diff_mat = np.abs(jax_mean - torch_mean)
    im = ax.imshow(diff_mat, cmap="Reds", aspect="equal")
    plt.colorbar(im, ax=ax, shrink=0.8)
    ax.set_xlabel("Sample")
    ax.set_ylabel("Sample")
    ax.set_title(f"|JAX - PyTorch| mean dist matrix\n(max={diff_mat.max():.2e})")

    # Panel 3: Histogram of element-wise differences
    ax = axes[2]
    diffs = jax_flat - torch_flat
    ax.hist(diffs, bins=100, color="#4CAF50", alpha=0.7, edgecolor="none")
    ax.axvline(x=0, color="k", linestyle="--", linewidth=0.8)
    ax.set_xlabel("JAX distance - PyTorch distance")
    ax.set_ylabel("Count")
    mean_abs = np.mean(np.abs(diffs))
    ax.set_title(f"Element-wise differences\n(mean |diff|={mean_abs:.2e})")
    ax.grid(True, alpha=0.3)

    fig.suptitle(
        "JAX vs PyTorch MRVI: Local Sample Distance Matrices (same weights, use_mean=True)",
        fontsize=13,
        y=1.02,
    )
    fig.tight_layout()
    fig.savefig(save_path, dpi=150, bbox_inches="tight")
    print(f"Saved distance comparison plot to {save_path}")
    plt.close(fig)


# ── Main ─────────────────────────────────────────────────────────────────────


if __name__ == "__main__":
    out_dir = os.path.join(os.path.dirname(__file__), "mrvi_equivalence_plots")
    os.makedirs(out_dir, exist_ok=True)

    print("=" * 72)
    print("1. PARAMETER UPDATE EQUIVALENCE")
    print("=" * 72)
    results = run_parameter_update_comparison(n_steps=15)
    plot_parameter_updates(results, os.path.join(out_dir, "parameter_updates.png"))

    print("\n" + "=" * 72)
    print("2. DISTANCE MATRIX + DE EQUIVALENCE")
    print("=" * 72)
    jax_de, torch_de, jax_dists, torch_dists, jax_reps, torch_reps = run_de_comparison()
    plot_de_comparison(
        jax_de, torch_de, jax_reps, torch_reps, os.path.join(out_dir, "de_comparison.png")
    )
    plot_distance_comparison(
        jax_dists, torch_dists, os.path.join(out_dir, "distance_comparison.png")
    )

    print("\n" + "=" * 72)
    print(f"All plots saved to {out_dir}/")
    print("=" * 72)
