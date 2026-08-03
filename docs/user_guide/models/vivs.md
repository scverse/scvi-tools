# VIVS

**VIVS** {cite:p}`BoyeauVIVS24` (Variational Inference for Variable Selection; Python class {class}`~scvi.external.VIVS`) identifies which genes in a count matrix `X` are conditionally dependent on an external response `Y` — protein expression, niche composition, or any other `obsm` feature — using a conditional randomization test (CRT).

VIVS fits two components: a generative VAE over `X` (reusing scvi-tools' standard `VAE` module, or an already-trained model passed via `x_model` — see "Compatible knockoff-sampler models" below), used only to sample conditional "knockoff" replacements of each gene; and an importance-score network predicting `Y` from `X`, whose per-cell negative log-likelihood is the CRT test statistic. For each gene (or, in the hierarchical variant, each cluster of correlated genes), the statistic is recomputed after substituting in a knockoff sample, yielding a calibrated, BH-corrected p-value for that gene's/cluster's conditional importance to `Y`.

The advantages of VIVS are:

-   Conditional (not marginal) dependence testing — controls for the rest of the transcriptome when assessing each gene's relevance to `Y`.
-   Calibrated false-discovery-rate control via the CRT + Benjamini-Hochberg correction.
-   A hierarchical variant (`get_hier_importance`) that tests at multiple resolutions of correlated gene clusters, improving power for co-regulated genes.
-   Can reuse an already-trained model (e.g. a niche-aware {class}`~scvi.external.SCVIVA` model) as the knockoff sampler — see the compatibility table below for which models qualify.

The limitations of VIVS include:

-   The knockoff sampler's quality (how well the generative VAE models `p(X_g | X_{-g})`) bounds the test's power.
-   Runtime scales with the number of genes tested; filtering to a smaller gene set (`select_genes`) is recommended above a few thousand genes.
-   Ported from VIVS's original JAX implementation; large-scale runtime is not guaranteed to match the original's `vmap`/`jit`-optimized performance (see `use_vmap` on `get_importance`/`get_hier_importance`).
-   `fastcluster` (for `get_gene_groupings`/`get_hier_importance`) and `plotnine` (for `plot_hier_importance`) are optional dependencies. `fastcluster` ships with `pip install scvi-tools[optional]`; `plotnine` has no dedicated lightweight extra — install it directly (`pip install plotnine`) or via `pip install scvi-tools[tutorials]`.

```{topic} Tutorials:

-   {doc}`/tutorials/notebooks/use_cases/VIVS_niche_gene_selection`
```

## Preliminaries

VIVS takes as input a raw-count gene expression matrix `X` and a response matrix `Y` (registered via `y_obsm_key` in {meth}`~scvi.external.VIVS.setup_anndata`). Training proceeds in two sequential phases: first the generative VAE over `X` is fit to convergence (or supplied pretrained via `x_model`), then it is frozen and the importance-score network is fit to predict `Y` from `X`. This order is required for CRT validity — the knockoff sampler must not be contaminated by information about `Y`.

## Compatible knockoff-sampler models

`VIVS(..., x_model=...)` requires `x_model.module` to be an instance of `scvi.module.VAE` (or a subclass), because the knockoff sampler calls `module.inference(x=x, batch_index=batch_index)` expecting a dict with `z`/`library`, then `module.generative(z=z, library=library, batch_index=batch_index)` expecting a dict with a `px` distribution supporting `.sample()`. Most non-`VAE` scvi-tools modules diverge from this signature (extra required arguments, or a different output-dict shape), so they fail either the `isinstance` check up front or the first knockoff-sampling call. Verified directly against module source (`src/scvi/module/`, `src/scvi/external/*/`):

| Model | `.module` class | Works as `x_model` today? | Notes |
|---|---|---|---|
| {class}`~scvi.model.SCVI` | `VAE` | ✅ Yes | Reference case; `x_module_kwargs` map onto `VAE.__init__`. |
| {class}`~scvi.external.SCVIVA` | `nicheVAE(VAE)` | ✅ Yes | Only overrides `__init__`/encoder-decoder wiring, not `inference`/`generative`. |
| {class}`~scvi.model.LinearSCVI` | `LDVAE(VAE)` | ✅ Yes | `inference`/`generative` fully inherited, unchanged. |
| {class}`~scvi.model.AUTOZI` | `AutoZIVAE(VAE)` | ✅ Yes | `generative` override only adds optional kwargs (dropout rescaling); `z`/`library`/`batch_index` call still works, still returns `px`. |
| {class}`~scvi.model.SCANVI` | `SCANVAE(SupervisedModuleClass, VAE)` | ✅ Yes | Same as `LinearSCVI` — no `inference`/`generative` override. |
| {class}`~scvi.external.GIMVI` | `JVAE(BaseModuleClass)` | ⚠️ Blocked only by the `isinstance` check | Not a `VAE` subclass, so `VIVS.__init__` currently rejects it — but its `inference()`/`generative()` happen to be duck-type compatible (`z`/`library` in, `px` out, `mode` defaults to `0`). Relaxing the `isinstance` gate to a duck-type/allowlist check would likely make this work with no other changes; untested. |
| {class}`~scvi.model.DestVI` | `MRDeconv(EmbeddingModuleMixin, BaseModuleClass)` | ❌ No | `generative(self, z, ind_x, library, batch_index)` — `ind_x` is a required positional argument VIVS never supplies; fails on the first knockoff call, not just the `isinstance` check. |
| {class}`~scvi.external.RESOLVI` | `RESOLVAE(PyroBaseModuleClass)` | ❌ No | Pyro-based, no `inference()`/`generative()` VAE API at all. |
| {class}`~scvi.model.TOTALVI` | `TOTALVAE(BaseMinifiedModeModuleClass)` | ❌ No | `generative` needs `library_gene` (not `library`) and a required `label`; output dict has no `px` key (`px_`/`py_` nested dicts instead). |
| {class}`~scvi.model.MULTIVI` | `MULTIVAE(BaseMinifiedModeModuleClass)` | ❌ No | `generative` needs a required `qz_m` positional arg; output dict has no `px` key (`px_scale`/`px_rate`/`px_dropout` separately). |
| {class}`~scvi.model.PEAKVI` | `PEAKVAE(BaseModuleClass)` | ❌ No | `generative` needs a required `qz_m` positional arg (not produced by VIVS's `_encode_for_knockoffs`), though it does return `px`. |
| {class}`~scvi.model.CondSCVI` | `VAEC(EmbeddingModuleMixin, BaseModuleClass)` | ❌ No | `generative` requires a cell-type-label `y` positional arg — semantically unrelated to (and would collide with) VIVS's own `y` (the response `Y`). |
| {class}`~scvi.model.AmortizedLDA` | `AmortizedLDAPyroModule(PyroBaseModuleClass)` | ❌ No | Pyro-based, no VAE-style `inference`/`generative`. |

Practical takeaway: today, the safest `x_model` choices are a plain {class}`~scvi.model.SCVI`, an {class}`~scvi.external.SCVIVA` model, or (with less real-world precedent) `LinearSCVI`/`AUTOZI`/`SCANVI`. Deconvolution (`DestVI`), Pyro-based (`RESOLVI`, `AmortizedLDA`), and paired/multimodal (`TOTALVI`, `MULTIVI`, `PEAKVI`, `CondSCVI`) models need either a model-specific adapter or a relaxed compatibility check before they can serve as VIVS's knockoff sampler.
