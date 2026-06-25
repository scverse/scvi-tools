# contrastiveVI

**contrastiveVI** {cite:p}`Weinberger23` (contrastive variational inference; Python class
{class}`~scvi.external.ContrastiveVI`) is a generative model for the contrastive analysis
of scRNA-seq count data. It is designed for experiments with a _background_ population,
such as control cells, and a _target_ population, such as perturbed or diseased cells.
After training, contrastiveVI can be used to extract target-specific latent structure,
estimate normalized expression, and perform differential expression.

Contrastive analysis requires a _target_ (e.g., treated cells) and a _background_
(e.g., control cells) dataset, and contrastiveVI is designed to isolate the variations
enriched in target cells from variations shared with background cells.

The advantages of contrastiveVI are:

-   It explicitly separates variation shared by target and background cells from variation
    enriched in the target cells.
-   It provides both background and salient latent representations, which can be used for
    visualization, clustering, and perturbation-specific downstream analyses.
-   It uses a count likelihood and supports many of the same data registration options as
    other scRNA-seq models in scvi-tools, including batches and observed covariates.

The limitations of contrastiveVI include:

-   It requires the user to define suitable background and target cell sets before training.
-   The salient latent space is interpretable relative to the supplied background population;
    if important target variation is also present in the background, contrastiveVI may not
    separate it as salient.
-   Like other nonlinear variational models, the individual latent dimensions are not
    guaranteed to be independently interpretable.

```{topic} Tutorials:

-   {doc}`/tutorials/notebooks/scrna/contrastiveVI_tutorial`
```

## Preliminaries

contrastiveVI takes as input a scRNA-seq count matrix $X$ with $N$ cells and $G$ genes.
The AnnData object is registered with {meth}`~scvi.external.ContrastiveVI.setup_anndata`;
the raw counts may be stored in `adata.X` or in a layer, and optional batch, label, and
covariate annotations may also be registered.

Training also requires two sets of cell indices:

-   `background_indices`, containing the background cells whose variation defines the
    shared baseline.
-   `target_indices`, containing the target cells in which salient variation should be
    isolated.

The training data splitter keeps separate background and target train, validation, and
test splits. Each contrastive minibatch contains one background batch and one target
batch, so the model learns both shared and target-enriched structure during each
optimization step.

## Model Overview

contrastiveVI uses two latent spaces:

-   $z_i$, the background latent variable, captures variation shared by background and
    target cells.
-   $s_i$, the salient latent variable, captures variation enriched in target cells.

The key modeling constraint is that background cells are decoded with $s_i = 0$, while
target cells are decoded with both $z_i$ and $s_i$. Both populations share the same
decoder. As a result, features that help reconstruct both background and target cells can
be represented in $z_i$, while variation needed specifically for target cells is encouraged
to enter $s_i$.

For a cell $i$, the decoder receives the concatenated latent representation
$[z_i, s_i]$, the library size, and registered batch information, and returns the
parameters of a zero-inflated negative binomial likelihood over genes. The model can use
the observed library size directly or infer a latent library size, controlled by
`use_observed_lib_size`.

## Inference

contrastiveVI uses amortized variational inference. Neural encoders approximate the
posterior distributions of the background latent variable $z_i$, salient latent variable
$s_i$, and, when needed, the library size. The training objective combines:

-   reconstruction terms for both background and target cells,
-   KL penalties for the background latent variables in both populations,
-   a KL penalty for the salient latent variables in target cells,
-   a library-size KL penalty when the library size is inferred, and
-   an optional Wasserstein penalty that discourages background-cell variation from
    leaking into the salient latent space.

In the implementation, the salient latent variables for background cells are explicitly
set to zero during inference before decoding.

## Tasks

Here we provide an overview of common tasks. Please see
{class}`~scvi.external.ContrastiveVI` for the full API reference.

### Salient and Background Representations

The default latent representation is the salient representation:

```
>>> adata.obsm["X_cvi_salient"] = model.get_latent_representation()
```

The shared background representation can be retrieved with:

```
>>> adata.obsm["X_cvi_background"] = model.get_latent_representation(
...     representation_kind="background"
... )
```

The salient representation is often the representation used for target-cell visualization
or clustering because it is trained to remove structure shared with background cells.

### Normalized Expression

{meth}`~scvi.external.ContrastiveVI.get_normalized_expression` returns decoded expression
from both latent spaces in a dictionary with `"background"` and `"salient"` entries. The
background output decodes the cell with the salient latent set to zero; the salient output
decodes the cell with the inferred salient latent included.

```
>>> exprs = model.get_normalized_expression()
>>> salient_expr = exprs["salient"]
>>> background_expr = exprs["background"]
```

For convenience, {meth}`~scvi.external.ContrastiveVI.get_salient_normalized_expression`
returns only the salient decoded expression.

### Differential Expression

Differential expression is available with
{meth}`~scvi.external.ContrastiveVI.differential_expression`. By default, the method uses
salient normalized expression, which focuses the comparison on target-enriched variation.
When `target_idx` is supplied, comparisons fully inside the target set use salient
expression, while comparisons outside the target set can fall back to the background
decoded expression.
