# Harreman

**Harreman** (`scvi.external.harreman`) is a toolkit for inferring metabolic exchanges and cell-cell communication in tissues using spatial transcriptomics data.

The advantages of Harreman are:

- Infers spatially-resolved metabolic gene programs using local autocorrelation
- Identifies cell-cell metabolic communication and ligand-receptor interactions using spatial proximity graphs
- Supports multiple spatial technologies (Visium, Slide-seq, and others)
- Scalable to large spatial datasets
- Supports both parametric and non-parametric significance testing

The limitations of Harreman include:

- Requires spatial coordinates to be available in `adata.obsm`
- Cell communication inference requires a ligand-receptor or metabolite transporter database

```{topic} Tutorials:

-   {doc}`/tutorials/notebooks/spatial/harreman_tutorial`
```

```{topic} External links:

- [Harreman documentation](https://harreman.readthedocs.io)
- [Harreman GitHub](https://github.com/YosefLab/Harreman)
```

## Overview

Harreman operates in three main steps:

1. **Spatial graph construction** ({func}`~scvi.external.harreman.tl.compute_knn_graph`): builds a spatial proximity graph from cell coordinates, supporting both k-nearest neighbors and radius-based neighborhoods, with optional Gaussian kernel weighting.

2. **Local autocorrelation** ({func}`~scvi.external.harreman.hs.compute_local_autocorrelation`): identifies spatially variable genes using the local autocorrelation statistic from the Hotspot algorithm (DeTomaso and Yosef, *Cell systems*, 2021), supporting DANB, Bernoulli, and normal count models.

3. **Cell communication** ({func}`~scvi.external.harreman.tl.compute_cell_communication`): infers spatially-resolved metabolic exchanges and ligand-receptor interactions between neighboring cells using HarremanDB and CellChatDB.

## Generative process

At the coarsest level, Harreman partitions the tissue into modules of different metabolic functions based on enzyme co-expression. At the following stage, Harreman formulates hypotheses about which metabolites are exchanged across the tissue or within each spatial zone. Moving to a finer resolution, Harreman can also infer which specific cell subsets participate in the exchange of distinct metabolic activities inside each zone.

For proteins composed of multiple subunits, Harreman computes either an algebraic or geometric mean of the expression values of the corresponding genes:

```{math}
:nowrap: true

\begin{align}
    X_{ai} &= \frac{\sum_{l \in S_l} X_{a_li}}{|S_l|}; \quad X_{bj} = \frac{\sum_{r \in S_r} X_{b_rj}}{|S_r|}
\end{align}
```

### Test statistic 1: Spatial autocorrelation

Spatially variable genes are identified using the following autocorrelation statistic:

```{math}
:nowrap: true

\begin{align}
    H_{a} &= \sum_{i}\sum_{j} w_{ij}X_{ai}X_{aj}
\end{align}
```

where $w_{ij}$ represents the communication strength between neighboring cells, computed using a Gaussian kernel:

```{math}
:nowrap: true

\begin{align}
    \hat{w}_{ij} &= e^{-d_{ij}^2/\sigma_{i}^2}
\end{align}
```

Significance is assessed by converting $H_a$ to a Z-score and adjusting p-values using the Benjamini-Hochberg procedure.

### Test statistic 2: Spatial co-localization

Pairwise spatial correlation between genes is computed as:

```{math}
:nowrap: true

\begin{align}
    H_{ab} &= \sum_{i}\sum_{j} w_{ij} \left(X_{ai}X_{bj} + X_{bi}X_{aj}\right)
\end{align}
```

This statistic is used to group genes into spatial modules and to identify cell-type-agnostic metabolic exchange events.

### Test statistic 3: Metabolite autocorrelation

Gene-pair results are integrated at the metabolite level:

```{math}
:nowrap: true

\begin{align}
    H_{m} &= \sum_{a,b \in m} H_{ab}
\end{align}
```

where $m$ is a metabolite exchanged by genes $a$ and $b$.

## Usage

```python
import scvi.external.harreman as harreman

# 1. Build spatial KNN graph
harreman.tl.compute_knn_graph(adata, compute_neighbors_on_key="spatial", n_neighbors=10)

# 2. Identify spatially variable genes
harreman.hs.compute_local_autocorrelation(adata, model="danb")

# 3. Compute pairwise local correlation
harreman.hs.compute_local_correlation(adata)

# 4. Infer cell-cell communication
harreman.tl.compute_cell_communication(adata)
```

## API

Please see {mod}`scvi.external.harreman` for the full API reference.
