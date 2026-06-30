# Tangram

:::{note}
This model is deprecated starting v1.5.
:::

**Tangram** {cite:p}`Biancalani21` (Python class {class}`~scvi.external.Tangram`) maps single-cell RNA-seq data to spatial data, permitting deconvolution of cell types in spatial data like Visium.

This is a reimplementation of Tangram, which can originally be found [here](https://github.com/broadinstitute/Tangram).

The advantages of Tangram are:

-   It maps single-cell transcriptomes onto spatial observations with a directly
    interpretable mapping matrix.
-   It can project cell annotations, such as cell types, from single-cell data to spatial
    data.
-   It can project gene expression from the single-cell reference into spatial
    coordinates.
-   The scvi-tools implementation supports Tangram's `"cells"` and `"constrained"`
    modes.

The limitations of Tangram include:

-   The scvi-tools model page is deprecated starting v1.5.
-   It requires matched genes between the single-cell and spatial modalities used for
    training.
-   Training is not mini-batched in the current implementation, so memory use depends on
    the number of single-cell observations and spatial observations.
-   Tangram is an optimization-based mapping model, not a generative model with posterior
    sampling.

## Overview

Tangram learns a matrix $M$ with shape ($n_{sc} \times n_{sp}$), in which each row sums to 1. Thus, this matrix can be viewed as a map from single cells to the spatial observations.

:::{note}
Starting scVI-Tools v1.5 this model is part of scVIVA-Tools, and no longer being maintained here.
:::

```{topic} Tutorials:

-   {doc}`/tutorials/notebooks/spatial/tangram_scvi_tools`
```

## Preliminaries

Tangram is registered with {meth}`~scvi.external.Tangram.setup_mudata`, not
`setup_anndata`. The input is a MuData object containing a single-cell modality and a
spatial modality. The two modalities used for training must contain the same genes in the
same order.

The `modalities` argument tells Tangram which MuData modality contains each registered
field. A typical setup registers:

-   `sc_layer`, the single-cell expression matrix used for mapping,
-   `sp_layer`, the spatial expression matrix used as the target, and
-   `density_prior_key`, an optional spatial observation column containing a density prior.

If a density prior is supplied, it must sum to 1. The tutorial computes a density prior
from estimated cell counts in spatial observations.

## Mapping Objective

The mapping matrix $M$ is parameterized by unconstrained trainable weights and converted
to a row-stochastic matrix with a softmax. The predicted spatial expression matrix is:

```{math}
:nowrap: true

\begin{align}
 \hat{X}_{sp} = M^\top X_{sc},
\end{align}
```

where $X_{sc}$ is the single-cell expression matrix and $\hat{X}_{sp}$ is the spatial
expression predicted from mapped single cells.

Tangram optimizes a loss that rewards agreement between measured spatial expression and
the predicted spatial expression. The default expression term uses gene-wise cosine
similarity, and an optional spatial-observation-wise cosine term can also be enabled. If a
density prior is registered, the loss includes a KL-divergence term between the predicted
spatial density and the supplied prior.

In constrained mode, Tangram also learns a cell filter and requires `target_count`. The
loss then includes a count term that encourages the selected number of cells to match
`target_count`, plus a filter regularizer.

## Tasks

Here we provide an overview of common tasks. Please see {class}`~scvi.external.Tangram`
for the full API reference.

### Mapping Matrix

After training, {meth}`~scvi.external.Tangram.get_mapper_matrix` returns the mapping
matrix with shape `(n_obs_sc, n_obs_sp)`:

```
>>> mapper = model.get_mapper_matrix()
```

Each row contains a probability distribution over spatial observations for one
single-cell observation.

### Projection of Cell Annotations

{meth}`~scvi.external.Tangram.project_cell_annotations` uses the mapping matrix to
project categorical single-cell labels, such as cell types, onto spatial observations:

```
>>> adata_sp.obsm["tangram_ct_pred"] = model.project_cell_annotations(
...     adata_sc, adata_sp, mapper, adata_sc.obs["cell_type"]
... )
```

The returned DataFrame has spatial observations as rows and label categories as columns.

### Projection of Genes

{meth}`~scvi.external.Tangram.project_genes` multiplies the mapping matrix by the
single-cell expression matrix and returns an AnnData object containing projected gene
expression in spatial coordinates.
