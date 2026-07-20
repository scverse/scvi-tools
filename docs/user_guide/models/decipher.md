# Decipher

**Decipher** {cite:p}`Nazaret24` (Python class {class}`~scvi.external.Decipher`) is a
probabilistic model for interpretable representation learning in single-cell RNA-seq
data. Decipher learns a low-dimensional latent representation for visualization and a
higher-dimensional latent representation for refined cell-state information, and it
provides utilities for imputed expression and trajectory-associated gene patterns.

The advantages of Decipher are:

-   It learns an interpretable latent space `v`, which is two-dimensional by default and
    can be visualized directly.
-   It also learns an intermediate latent space `z`, which can capture more detailed
    cell-state structure than the visualization space.
-   It includes helper methods for imputed gene expression, Decipher time, and gene
    expression patterns along a trajectory.

The limitations of Decipher include:

-   The current scvi-tools implementation registers only the count matrix, with an
    optional raw-count layer, and does not yet expose condition covariates through
    {meth}`~scvi.external.Decipher.setup_anndata`.
-   The tutorial currently demonstrates the basic Decipher model; a fuller implementation
    aligned with the original method is still under development.
-   As with other latent-variable models, trajectory interpretation depends on the quality
    of the learned latent representation and the user-specified trajectory.

```{topic} Tutorials:

-   {doc}`/tutorials/notebooks/scrna/decipher_tutorial`
```

## Preliminaries

Decipher takes as input a scRNA-seq count matrix $X$ with $N$ cells and $G$ genes. The
AnnData object is registered with {meth}`~scvi.external.Decipher.setup_anndata`; by
default, the model uses `adata.X`, but a count layer can be provided with the `layer`
argument.

The model does not require batch annotations or other covariates in the current
implementation. After registration, a Decipher model can be created and trained with the
usual scvi-tools pattern:

```
>>> Decipher.setup_anndata(adata, layer="counts")
>>> model = Decipher(adata)
>>> model.train()
```

## Model Overview

Decipher uses two latent representations:

-   $v_i$, a low-dimensional interpretable representation of cell $i$.
-   $z_i$, an intermediate latent representation that links $v_i$ to gene expression.

The generative model starts from a standard normal prior on $v_i$. A decoder maps
$v_i$ to the parameters of a normal distribution over $z_i$, and a second decoder maps
$z_i$ to gene proportions. These gene proportions are multiplied by the observed library
size of the cell and used as the mean of a negative binomial observation model with
learned gene-specific inverse dispersion.

The default dimensions are `dim_v=2` and `dim_z=10`. The two-dimensional `v`
representation is intended for direct visualization, while `z` is intended to retain
more information for reconstruction and downstream cell-state analyses.

## Inference

Decipher is implemented as a Pyro module and trained with stochastic variational
inference. The guide uses encoder networks in the reverse direction: log-transformed
counts are encoded to $z_i$, and the concatenation of $z_i$ and the log-transformed
counts is encoded to $v_i$. The `beta` parameter scales the KL term for $v_i`, which
controls the regularization strength of the interpretable latent space.

During training, Decipher uses validation predictive log likelihood as the default
early-stopping monitor when early stopping is enabled.

## Tasks

Here we provide an overview of common tasks. Please see {class}`~scvi.external.Decipher`
for the full API reference.

### Latent Representation

The default latent representation is $v_i$:

```
>>> adata.obsm["X_decipher_v"] = model.get_latent_representation()
```

The intermediate representation $z_i$ can be returned with `give_z=True`:

```
>>> adata.obsm["X_decipher_z"] = model.get_latent_representation(give_z=True)
```

### Imputed Expression

{meth}`~scvi.external.Decipher.compute_imputed_gene_expression` decodes cells through the
model and returns imputed gene expression on the scale of each cell's observed library
size. When `compute_covariances=True`, the method also returns covariances between the
imputed expression and the stored Decipher `v` and `z` representations.

### Trajectories and Gene Patterns

The Decipher utilities include a `Trajectory` object for representing paths through the
Decipher latent space. Given a trajectory and cluster assignments,
{meth}`~scvi.external.Decipher.compute_decipher_time` estimates a Decipher time for cells
on the trajectory with K-nearest-neighbor regression. The
{meth}`~scvi.external.Decipher.compute_gene_patterns` method then decodes points along
the trajectory to summarize gene expression patterns with means and quantiles.
