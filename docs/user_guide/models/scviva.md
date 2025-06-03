# scVIVA

**scVIVA** (Python class {class}`~scvi.external.SCVIVA`) is a generative model of single-cell resolved spatial
transcriptomics that can subsequently be used for many common downstream tasks.

The advantages of scVIVA are:

-   Provides a probabilistic low-dimensional representation of the state of each cell that is corrected for batch effects
    and captures its gene expression profile and its environment.
-   Enables differential expression analysis across niches while accounting for wrong assignment of molecules to cells.
-   Scalable to very large datasets (>1 million cells).

The limitations of scVIVA include:

-   Effectively requires a GPU for fast inference.
-   Latent space is not interpretable, unlike that of a linear method.
-   Assumes single cells are observed and does not work with low resolution ST like Visium or Slide-Seq.

```{topic} Tutorials:

-   {doc}`/tutorials/notebooks/spatial/scVIVA_tutorial`
```

## Preliminaries

scVIVA takes as input spatially-resolved scRNA data. In addition to the gene expression matrix ${X}$ with $N$ cells and $G$ genes,
it requires for each cell $n$:
- the spatial coordinates of the cell $y_n$
- the cell type assignment (possibly coarse) $c_n \in \{1, ..., T\}$
- the batch assignment $s_n$.


As preprocessing, we take the $K$ nearest neighbors of a cell to define its niche using the Euclidean distance in physical space.
We characterize the niche by its cell-type composition and gene expression. We denote by ${\alpha_n}$ the $T$-dimensional vector of cell type
proportions among the $K$ nearest neighbors of the cell $n$. Its values are in the probability simplex.
The niche gene expression is defined as the average expression of each cell type present in the niche.
In practice, we leverage gene expression embeddings (PCA, scVI or similar) and characterize a cell type expression profile as the local average
embedding of cells of the same type. The average embeddings are stored in the matrix ${\eta_n} \in \mathbb{R}^{T \times D}$, where $D$ is the embedding dimension.
## Descriptive model

We propose a latent variable model aiming to capture both gene expression heterogeneity and spatial variation resulting from the micro-environment.
We assume these two sources of variability are both captured by a $P$-dimensional latent variable ($P \ll G$):

```{math}
:nowrap: true
\begin{align}
    z_n \sim \mathbf{MixtureOfGaussians}(\mu_1, ..., \mu_M; \Sigma_1, ..,\Sigma_M; \pi_1, ...,\pi_M)
\end{align}
```

We assume that the observed counts for cell $n$ and gene $g$, $x_{ng}$, are generated from the following process:

```{math}
:nowrap: true
\begin{align}
 \rho _n &= f_{w}\left( z_n, s_n \right) \\
 x_{ng} &\sim \mathbf{NegativeBinomial}(\ell_n \rho_n, \theta_g),
 \end{align}
```
where $\rho_n$ is the normalized gene expression, $\ell_n$ is the library size of cell $n$ and $\theta_g$ is the dispersion parameter for gene $g$. See {doc}`/user_guide/models/scvi` for more details.
The cell-type proportions of the cell's $K$ nearest neighbors are obtained as

```{math}
:nowrap: true
\begin{align}
    \alpha_n &\sim \mathbf{Dirichlet}\left( f_{\omega}(z_n) \right),
\end{align}
```

Last, we assume that the neighboring cells' average expression profiles are obtained as

```{math}
:nowrap: true
\begin{equation}
\eta_{nt} \sim
\begin{cases}
\mathcal{N} \left(f_{\nu}^{t}(z_n) \right), & \text{if } \alpha_{t} > 0 \\
0, & \text{otherwise}
\end{cases}
\end{equation}
```

where $t=1,...,T$. $w$, $\omega$ and $\nu$ are neural network parameters.


## Inference

We want to maximize the evidence of the data, which can be decomposed as:

```{math}
:nowrap: true
\begin{align}
    \log p \left( \alpha, x, \eta \mid s \right) = \log p \left(x \mid s \right) + \log p \left( \alpha, \eta \mid x, s \right).
\end{align}
```

scVIVA uses variational inference, specifically auto-encoding variational Bayes
(see {doc}`/user_guide/background/variational_inference`) to learn both the model parameters
(the neural network parameters, dispersion parameters, etc.) and an approximate posterior distribution.

## Tasks

Here we provide an overview of some of the tasks that scVIVA can perform. Please see {class}`scvi.external.SCVIVA`
for the full API reference.

### Dimensionality reduction

For dimensionality reduction, the mean of the approximate posterior $q_\phi(z \mid x)$ is returned by default.
This is achieved using the method:

```
>>> adata.obsm["X_scVIVA"] = model.get_latent_representation()
```

$\phi$ is a set of parameters corresponding to inference neural networks (encoders).
Users may also return samples from this distribution, as opposed to the mean, by passing the argument `give_mean=False`.

### Estimation of normalized expression

In {meth}`~scvi.external.SCVIVA.get_normalized_expression` scVIVA returns the expected true expression value of $\rho_n$ under the approximate posterior. For one cell $n$, this can be written as:

```{math}
:nowrap: true

\begin{align}
   \mathbb{E}_{q_\phi(z_n \mid x_n)}\left[f_{w}\left(z_{n}, s_n \right) \right]
\end{align}
```

### Differential Expression (DE)

Differential expression analysis is achieved with {meth}`~scvi.external.SCVIVA.differential_expression`. \
We leverage the lvm-DE method (see {doc}`/user_guide/background/differential_expression`) and adapt it to spatial data by taking into account cell neighborhood expression in a bid to discard false positives due to contamination. \
Considering two groups of cells $\textit{C1}$ and $\textit{C2}$ corresponding to different spatial contexts (for instance, astrocytes in two brain regions), the goal is to determine which genes have different expression levels between the two groups. When setting `niche_mode="true"`, we compute the group spatial neighborhoods $\textit{N1}$ and $\textit{N2}$, which are the spatial nearest neighbors of a different type than the cells in $\textit{C1}$, and $\textit{C2}$ respectively.


To determine the upregulated genes of $\textit{C1 vs C2}$, we compute DE between $\{\textit{C1, C2}\}$, $\{\textit{N1, C2}\}$ and $\{\textit{C1, N1}\}$: using lvm-DE, we test differences in expression levels $\rho_{n}$ to compute Log-Fold Changes (LFC). \
The significantly upregulated genes for $\textit{C1, N1}$ define a set of local cell type markers, denoted $\mathcal{S}_1$. Conversely, if a gene is both higher expressed in $\textit{N1}$ compared to $\textit{C1}$ and $\textit{C1}$ compared to $\textit{C2}$, it is likely that the increased expression in $\textit{C1}$ is spurious.
We argue that the probability of a gene being a $\textit{local marker}$ could be a relevant score to filter spurious genes. To compute this score, we considered the upregulation of a gene in one group relative to the upregulation in its neighborhood: a local marker $g$ should verify

```{math}
:nowrap: true
\begin{align}
    \mathit{LFC^{~g}_{C1~vs~C2}} > \mathit{LFC^{~g}_{N1~vs~C2}},
\end{align}
```

which means that the signal comes from cells in $\textit{C1}$ rather than their neighbors $\textit{N1}$. \
We select genes for which $\mathit{LFC_{C1~vs~C2}} > 0$ and use the genes $\mathcal{S}_1$ as truely differentially expressed. We also define $\mathcal{N}_1 = \{g|\mathit{LFC^{~g}_{C1~vs~C2}} > 0,~g \notin \mathcal{S}_1 \}$. \
We train a Gaussian process classifier on $\mathbf{X} = [LFC_{C1~vs~C2}~,~LFC_{N1~vs~C2}]$  to classify between the $\textit{local markers}$ $\mathcal{S}_1$ and the $\textit{neighborhood genes}$ $\mathcal{N}_1$. Once fitted, the classifier returns a local marker probability $p_g=\mathit{p}(g \in \mathcal{S}_1 | \mathbf{X})$ for each gene $g$, that we can compare to a given threshold $\tau$ to filter the neighborhood genes.
