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
    adata = scvi.data.synthetic_iid()
    n_obs = adata.n_obs
    adata.obs["sample"] = np.random.choice(5, n_obs).astype(str)
    meta_per_sample = {str(i): np.random.randint(0, 2) for i in range(5)}
    adata.obs["meta1"] = adata.obs["sample"].map(meta_per_sample).astype(int)
    adata.obs["meta1_cat"] = "CAT_" + adata.obs["meta1"].astype(str)
    adata.obs["meta1_cat"] = adata.obs["meta1_cat"].astype("category")

    # ── Train JAX model ──
    print("\nTraining JAX model...")
    MRVI.setup_anndata(adata, sample_key="sample", batch_key="batch", backend="jax")
    jax_model = MRVI(adata)
    jax_model.train(max_epochs=50, batch_size=256, train_size=0.9)
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
    # DE uses stochastic MC sampling (use_mean=False internally), so results will
    # differ slightly due to different RNG systems. We use high mc_samples to minimize
    # variance, and also store_lfc=True with more training for meaningful LFCs.
    de_kwargs = {
        "sample_cov_keys": ["meta1_cat"],
        "store_lfc": True,
        "mc_samples": 250,
        "batch_size": 400,
        "use_vmap": False,
    }
    print("Running JAX DE (mc_samples=250)...")
    jax_de = jax_model.differential_expression(**de_kwargs)

    print("Running PyTorch DE (mc_samples=250)...")
    torch_de = torch_model.differential_expression(**de_kwargs)

    # Also run a second PyTorch DE with different seed to establish MC noise baseline
    print("Running PyTorch DE (mc_samples=250, different seed for baseline)...")
    torch.manual_seed(99)
    torch_de_2 = torch_model.differential_expression(**de_kwargs)

    return jax_de, torch_de, torch_de_2, jax_dists, torch_dists


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


def plot_de_comparison(jax_de, torch_de, torch_de_2, save_path):
    """Plot DE comparison between JAX and PyTorch with identical model weights.

    Both models share the same trained weights but use different RNG systems for
    MC sampling. Also shows a PyTorch-vs-PyTorch baseline (different seeds) to
    demonstrate that the JAX-vs-PyTorch spread matches the inherent MC noise.
    """
    from scipy.stats import spearmanr

    fig, axes = plt.subplots(1, 4, figsize=(20, 4.5))

    # Compute PyTorch baseline correlations
    t1_beta = torch_de["beta"].values.flatten()
    t2_beta = torch_de_2["beta"].values.flatten()
    baseline_beta_corr = np.corrcoef(t1_beta, t2_beta)[0, 1]

    # ── Panel 1: Beta coefficients scatter ──
    ax = axes[0]
    jax_beta = jax_de["beta"].values.flatten()
    torch_beta = torch_de["beta"].values.flatten()
    ax.scatter(jax_beta, torch_beta, s=8, alpha=0.5, color="#2196F3", edgecolors="none")
    lim = max(np.abs(jax_beta).max(), np.abs(torch_beta).max()) * 1.1
    ax.plot([-lim, lim], [-lim, lim], "k--", linewidth=0.8, alpha=0.5)
    ax.set_xlabel("JAX beta")
    ax.set_ylabel("PyTorch beta")
    corr = np.corrcoef(jax_beta, torch_beta)[0, 1]
    ax.set_title(
        f"Beta coefficients\n"
        f"JAX vs Torch: r={corr:.3f}\n"
        f"Torch vs Torch baseline: r={baseline_beta_corr:.3f}"
    )
    ax.set_xlim(-lim, lim)
    ax.set_ylim(-lim, lim)
    ax.set_aspect("equal")
    ax.grid(True, alpha=0.3)

    # ── Panel 2: Effect sizes scatter ──
    ax = axes[1]
    jax_es = jax_de["effect_size"].values.flatten()
    torch_es = torch_de["effect_size"].values.flatten()
    t2_es = torch_de_2["effect_size"].values.flatten()
    baseline_es_corr = np.corrcoef(torch_es, t2_es)[0, 1]
    ax.scatter(jax_es, torch_es, s=8, alpha=0.5, color="#4CAF50", edgecolors="none")
    lim_es = max(jax_es.max(), torch_es.max()) * 1.1
    ax.plot([0, lim_es], [0, lim_es], "k--", linewidth=0.8, alpha=0.5)
    ax.set_xlabel("JAX effect size")
    ax.set_ylabel("PyTorch effect size")
    corr_es = np.corrcoef(jax_es, torch_es)[0, 1]
    ax.set_title(
        f"Effect sizes\n"
        f"JAX vs Torch: r={corr_es:.3f}\n"
        f"Torch vs Torch baseline: r={baseline_es_corr:.3f}"
    )
    ax.grid(True, alpha=0.3)

    # ── Panel 3: LFC scatter ──
    ax = axes[2]
    if "lfc" in jax_de and "lfc" in torch_de:
        jax_lfc = jax_de["lfc"].values.flatten()
        torch_lfc = torch_de["lfc"].values.flatten()
        t2_lfc = torch_de_2["lfc"].values.flatten()
        mask = np.isfinite(jax_lfc) & np.isfinite(torch_lfc)
        mask_b = np.isfinite(torch_lfc) & np.isfinite(t2_lfc)
        baseline_lfc_corr = np.corrcoef(torch_lfc[mask_b], t2_lfc[mask_b])[0, 1]
        jax_lfc_m, torch_lfc_m = jax_lfc[mask], torch_lfc[mask]
        ax.scatter(jax_lfc_m, torch_lfc_m, s=4, alpha=0.3, color="#FF5722", edgecolors="none")
        lim_lfc = max(np.abs(jax_lfc_m).max(), np.abs(torch_lfc_m).max()) * 1.1
        ax.plot([-lim_lfc, lim_lfc], [-lim_lfc, lim_lfc], "k--", linewidth=0.8, alpha=0.5)
        ax.set_xlabel("JAX LFC")
        ax.set_ylabel("PyTorch LFC")
        corr_lfc = np.corrcoef(jax_lfc_m, torch_lfc_m)[0, 1]
        ax.set_title(
            f"Log fold changes\n"
            f"JAX vs Torch: r={corr_lfc:.3f}\n"
            f"Torch vs Torch baseline: r={baseline_lfc_corr:.3f}"
        )
        ax.set_xlim(-lim_lfc, lim_lfc)
        ax.set_ylim(-lim_lfc, lim_lfc)
        ax.set_aspect("equal")
    else:
        ax.text(
            0.5,
            0.5,
            "LFC not available",
            ha="center",
            va="center",
            transform=ax.transAxes,
            fontsize=12,
        )
        ax.set_title("Log fold changes")
    ax.grid(True, alpha=0.3)

    # ── Panel 4: P-value comparison ──
    ax = axes[3]
    jax_pv = jax_de["pvalue"].values.flatten()
    torch_pv = torch_de["pvalue"].values.flatten()
    t2_pv = torch_de_2["pvalue"].values.flatten()
    baseline_pv_rho = spearmanr(torch_pv, t2_pv).statistic
    ax.scatter(jax_pv, torch_pv, s=8, alpha=0.5, color="#9C27B0", edgecolors="none")
    ax.plot([0, 1], [0, 1], "k--", linewidth=0.8, alpha=0.5)
    ax.set_xlabel("JAX p-value")
    ax.set_ylabel("PyTorch p-value")
    rho_pv = spearmanr(jax_pv, torch_pv).statistic
    ax.set_title(
        f"P-values (Spearman)\n"
        f"JAX vs Torch: rho={rho_pv:.3f}\n"
        f"Torch vs Torch baseline: rho={baseline_pv_rho:.3f}"
    )
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.set_aspect("equal")
    ax.grid(True, alpha=0.3)

    fig.suptitle(
        "JAX vs PyTorch MRVI: Differential Expression (same weights)",
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
    jax_de, torch_de, torch_de_2, jax_dists, torch_dists = run_de_comparison()
    plot_de_comparison(jax_de, torch_de, torch_de_2, os.path.join(out_dir, "de_comparison.png"))
    plot_distance_comparison(
        jax_dists, torch_dists, os.path.join(out_dir, "distance_comparison.png")
    )

    print("\n" + "=" * 72)
    print(f"All plots saved to {out_dir}/")
    print("=" * 72)
