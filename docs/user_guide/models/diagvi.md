# DiagVI

**DiagVI** (Diagonal Multi-Modal Integration Variational Inference; Python class {class}`~scvi.external.DIAGVI`) is a deep generative model for integrating unpaired multi-modal single-cell data using prior biological knowledge encoded as a guidance graph.

DiagVI is inspired by the GLUE architecture , which uses modality-specific variational autoencoders (VAEs) to project heterogeneous data types into a shared latent space. In contrast to GLUE’s adversarial alignment strategy, DiagVI aligns modalities using Unbalanced Optimal Transport (UOT) via the Sinkhorn divergence, explicitly accounting for differences in cell-type composition across modalities.

This design makes DiagVI particularly well suited for settings with modality-specific or rare cell populations, such as spatial proteomics vs. scRNA-seq.

The advantages of DiagVI are:

-   Flexible two-modality integration of various data types (e.g., scRNA-seq, spatial transcriptomics, spatial proteomics).
-   Modality specific likelihoods tailored to the properties of the data.
-   Full feature utilization: all features (not only overlapping ones) contribute to model training via the guidance graph.
- Biologically informed alignment using prior feature correspondences.
-   Optional Gaussian mixture prior and semi-supervised learning with cell-type labels.

The limitations of DiagVI include:

-   Supports two modalities only.
-   Requires prior information on cross-modal feature correspondences (explicitly or implicitly).

```{topic} Tutorials:
 Tutorial placeholder
```

## Preliminaries

DiagVI takes as input two unpaired modalities with independent feature sets $\mathcal{V}_1, \mathcal{V}_2$ and observations $N_1, N_2$. 

### Supported modalities
Currently supported modalities include:

- scRNA-seq
- Spatial transcriptomics
- Spatial proteomics
- Other count or continuous measurements

Each modality may optionally include a batch covariate.

### Guidance Graph
Feature correspondences are encoded in the guidance graph $\mathcal{G} = (\mathcal{V}, \mathcal{E})$, where $\mathcal{V} = \mathcal{V}_1 \cup \mathcal{V}_2$.
Each edge $\mathcal{E} = \{(i,j)\mid i,j \in \mathcal{V}\}$ is associated with
- a weight $w_{ij} \in (0, 1]$ reflecting the confidence of the link
- a sign $s_{ij} \in [-1, 1]$ specifying whether the interaction is associative ($s_{ij} = 1$) or repressive ($s_{ij} = -1$)

The guidance graph regularizes feature embeddings to ensure biologically meaningful alignment between modalities. There are three ways to specify it:

- Automatic construction (default): \
If no graph is provided, DiagVI constructs a graph from overlapping feature names across modalities (e.g., shared gene symbols).

- Custom mapping via DataFrame: \
Useful when feature naming conventions differ (e.g., genes vs. proteins in CITE-seq). The DataFrame specifies feature pairs and optional weights/signs.

- Explicit graph specification: \
For full control, pass a `torch_geometric.data.Data` object directly.

## Descriptive model
<span style="color:red"> TODO: plate model </span>

DiagVI assumes that observations from each modality are generated from a shared 
$m$-dimensional latent space. 
For cell $n$ and feature $i$, the observed data $\mathbf{x}_{ni}$ is modelled as a function of the cell latent variable $\mathbf{z}_n \in \mathbb{R}^{m}$ and the feature latent variable  $\mathbf{v}_i \in \mathbb{R}^{m}$.

The prior on the feature (graph) latent is set to a standard normal distribution.
$$
    p(\boldsymbol{V}) \sim \prod_{i\in \mathcal{V}} N(\mathbf{v}_i;0,\boldsymbol{I_m})

$$

The prior on the cell latent variable can be either a standard Gaussian (default) or a Gaussian mixture with $L$ components.
$$
\mathbf{z}_{n} \sim
\begin{cases}
\mathcal{N}(0, 1), \\
\sum_{l=1}^L \pi_{l} \, \mathcal{N}(\mu_{l}, \sigma_{l}^2)
\end{cases}
$$
If cell-type labels are provided, a Gaussian mixture prior is used with one component per cell type. A modality-specific classifier predicting labels from 
$\mathbf{z}$ is trained jointly.

## Generative process
The graph decoder encourages embeddings $v_i$ and $v_j$ of linked features to be close in latent space, with strength modulated by edge weights and signs. Unrelated features are pushed apart.

Given cell latent $\mathbf{z}_n$ and feature embeddings $\mathbf{V}$, DiagVI generates the denoised, normalized data
$$
\boldsymbol{\rho}_n = \mathrm{softmax}\left( \boldsymbol{\alpha}_n \odot \left(  \mathbf{z}_n\mathbf{V}^\top \right) + \boldsymbol{\beta}_n \right)
$$
where $\boldsymbol{\alpha}_n \in \mathbb{R}_+^{\mid \mathcal{V} \mid}$ and $\boldsymbol{\beta}_n \in \mathbb{R}_+^{\mid \mathcal{V} \mid}$ are feature specific scaling and bias parameters, respectively. When a batch covariate is provided, batch-specific versions of these parameters are learned.

Observed library sizes $l_n$ are used to reconstruct raw counts.

## Likelihood models

Depending on the modality, different likelihood functions can be used to reconstruct the input data and the respective parameters are learned, e.g., feature-specific dispersion parameter $\boldsymbol{\theta} \in \mathbb{R}_+^{\mid \mathcal{V} \mid}$

$$
 x_{ni} \sim \text{NegativeBinomial} \left(l_n\rho_{ni}, \theta_i \right)
$$

---
DiagVI supports the following likelihood functions:

<span style="color:red">TODO: more description (maybe add a data preprocessing section?) </span>
```{eval-rst}
.. list-table::
   :widths: 20 40 40
   :header-rows: 1

   * - Likelihood
     - Distribution
     - modality type 
   * - ``nb``
     - Negative Binomial
     - scRNA-seq
   * - ``zinb``
     - Zero-Inflated Negative Binomial
     - scRNA-seq
   * - ``nbmixture``
     - Negative Binomial Mixture
     - Protein data (background/foreground)
   * - ``normal``
     - Normal
     - Protein data
   * - ``lognormal``
     - Log-Normal
     - Positive continuous data
   * - ``log1pnormal``
     - Log1p-Normal
     - Non-negative continuous data
   * - ``ziln``
     - Zero-Inflated Log-Normal
     - Sparse positive continuous data
   * - ``zig``
     - Zero-Inflated Gamma
     - Sparse positive continuous data
```

### Latent variables
<span style="color:red">needs some overthinking</span>
```{eval-rst}
.. list-table::
   :widths: 20 70 15
   :header-rows: 1

   * - Latent variable
     - Description
     - Code variable
   * - :math:`z_n \in \mathbb{R}^m`
     - Low-dimensional cell representation capturing biological state
     - ``z``
   * - :math:`v_i \in \mathbb{R}^{d}`
     - Low-dimensional feature embedding learned via the graph encoder
     - ``v``
   * - :math:`\rho_n \in \Delta^{G-1)}` or :math:`\mathbb{R}^G_+`
     - Normalized/denoised expression (sums to 1 for count data)
     - ``px_scale``
   * - :math:`\ell_n \in (0, \infty)`
     - Observed library size (total counts per cell)
     - ``library``
   * - :math:`\theta_g \in (0, \infty)`
     - Inverse dispersion parameter (per gene/batch)
     - ``px_r``
```
Depending on the likelihood and whether a GMM prior is used, these latent variables can vary.

## Inference

DiagVI uses variational inference (see {doc}`/user_guide/background/variational_inference`) to learn model parameters and approximate posterior distributions.

### Unbalanced Optimal Transport alignment
To align modality-specific latent spaces, DiagVI minimizes the de-biased Sinkhorn divergence between latent distributions using the GeomLoss library.

Key hyperparameters:
- `blur`: entropic regularization controlling transport smoothness,
- `reach`: marginal relaxation penalizing deviations from mass conservation.

This formulation naturally handles unequal cell numbers and modality-specific populations.

### Training objective

The total loss is a weighted sum of:

- Graph reconstruction loss (`lam_graph`),
- Data reconstruction loss (`lam_data`),
- KL divergence (`lam_kl`),
- UOT alignment loss (`lam_sinkhorn`),
- Optional classification loss (`lam_class`).

## Practical guidance
### Hyperparameter Selection
<span style="color:red"> TODO: Provide some recommendations for parameter selection. </span>


## Tasks

Here we provide an overview of the main tasks that DiagVI can perform. See {class}`~scvi.external.DIAGVI` for the full API reference.

### Dimensionality reduction
```python
>>> latents = model.get_latent_representation()
>>> adata_rna.obsm["X_diagvi"] = latents["rna"]
>>> adata_protein.obsm["X_diagvi"] = latents["protein"]
```
Aligned representations can be used jointly for downstream analysis:

```python
>>> import scanpy as sc
>>> # Use RNA latent for clustering
>>> sc.pp.neighbors(adata_rna, use_rep="X_diagvi")
>>> sc.tl.umap(adata_rna)
>>> sc.tl.leiden(adata_rna)
```

### Cross-modal feature imputation

DiagVI can impute features from one modality to another using {func}`~scvi.external.DIAGVI.get_imputed_values`.

```python
>>> # Impute protein expression from RNA data
>>> imputed_protein = model.get_imputed_values(
...     source_name="rna",
...     source_adata=adata_rna,
... )
>>> adata_rna.obsm["imputed_protein"] = imputed_protein
```

You can also specify target batch and library size for counterfactual predictions:

```python
>>> # Impute with specific target batch
>>> imputed = model.get_imputed_values(
...     source_name="rna",
...     target_batch="batch_1",
...     target_libsize=10000,
... )
```

## References

[^ref1]:
    Zhi-Jie Cao, Ge Gao (2022),
    _Multi-omics single-cell data integration and regulatory inference with graph-linked embedding_,
    [Nature Biotechnology](https://www.nature.com/articles/s41587-022-01284-4).
