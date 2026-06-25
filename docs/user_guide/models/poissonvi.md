# PoissonVI

**PoissonVI** {cite:p}`Martens2023` (Python class {class}`~scvi.external.POISSONVI`) is
a variational model for quantitative scATAC-seq fragment counts. It follows the
peakVI-style single-cell ATAC workflow, but uses a Poisson observation model to analyze
fragment counts rather than binarized accessibility.

The advantages of PoissonVI are:

-   It models quantitative fragment-count information, which can retain signal that is
    discarded when scATAC-seq matrices are binarized.
-   It exposes familiar peakVI-style tasks, including latent representation learning,
    normalized accessibility, region factors, and differential accessibility.
-   It supports registered batch and covariate annotations through
    {meth}`~scvi.external.POISSONVI.setup_anndata`.

The limitations of PoissonVI include:

-   It expects fragment-count-like input. Read-count matrices with strong even-count
    artifacts should be converted to fragment counts before training.
-   It is designed for chromatin accessibility regions, so
    {meth}`~scvi.external.POISSONVI.differential_expression` is not implemented; use
    {meth}`~scvi.external.POISSONVI.differential_accessibility` instead.
-   The nonlinear latent space is useful for integration and clustering, but individual
    latent dimensions are not directly interpretable.

```{topic} Tutorials:

-   {doc}`/tutorials/notebooks/atac/PoissonVI`
```

## Preliminaries

PoissonVI takes as input a cell-by-region matrix $X$ with $N$ cells and $M$ genomic
regions. The values should be quantitative fragment counts. For 10x scATAC-seq data, the
tutorial demonstrates converting read counts to approximate fragment counts with
{func}`~scvi.data.reads_to_fragments` and storing them in a layer before calling
{meth}`~scvi.external.POISSONVI.setup_anndata`.

The registered AnnData object may include a batch key, labels, size factors, and other
categorical or continuous covariates:

```
>>> POISSONVI.setup_anndata(adata, layer="fragments", batch_key="batch")
>>> model = POISSONVI(adata)
>>> model.train()
```

The setup method checks that the registered matrix has fragment-count-like behavior.

## Model Overview

PoissonVI reuses the variational autoencoder machinery used by other scvi-tools models,
with defaults chosen to mirror peakVI where appropriate. A cell-level latent variable
$z_i$ captures biological variation in chromatin accessibility. The decoder maps this
latent representation, along with registered covariates, to region-specific Poisson
rates.

In the implementation, PoissonVI uses {class}`~scvi.module.VAE` with
`gene_likelihood="poisson"`, gene-level dispersion settings that are unused by the
Poisson likelihood, layer normalization in both encoder and decoder, and LeakyReLU
activation functions.

## Inference

PoissonVI uses amortized variational inference. The encoder approximates the posterior
distribution over the latent accessibility representation, and the decoder is optimized
to reconstruct the registered fragment-count matrix under a Poisson likelihood.

The model also learns region-specific decoder factors. By default,
{meth}`~scvi.external.POISSONVI.get_normalized_accessibility` temporarily removes the
region-factor bias when returning normalized accessibility, which emphasizes
cell-specific accessibility variation. Passing `normalize_regions=True` reintroduces the
region factors and produces estimates closer to the observed input scale.

## Tasks

Here we provide an overview of common tasks. Please see {class}`~scvi.external.POISSONVI`
for the full API reference.

### Dimensionality Reduction

The latent representation can be stored in `adata.obsm` and used with Scanpy for
neighbors, clustering, and visualization:

```
>>> adata.obsm["X_poissonvi"] = model.get_latent_representation()
```

### Normalized Accessibility

{meth}`~scvi.external.POISSONVI.get_normalized_accessibility` returns normalized
accessibility estimates for cells and regions:

```
>>> accessibility = model.get_normalized_accessibility()
```

The method supports posterior sampling, region subsetting, `transform_batch`, and
importance weights. The learned region factors can be retrieved with
{meth}`~scvi.external.POISSONVI.get_region_factors`.

### Differential Accessibility

Differential accessibility is performed with
{meth}`~scvi.external.POISSONVI.differential_accessibility`. The method supports both
`"change"` and `"vanilla"` modes, one-sided or two-sided tests, and optional batch
correction.
