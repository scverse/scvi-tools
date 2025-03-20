# NicheVI

**NicheVI** (Python class {class}`~scvi.external.NICHEVI`) is a generative model of single-cell resolved spatial
transcriptomics that can subsequently be used for many common downstream tasks.

The advantages of NicheVI are:

-   Provides a probabilistic low-dimensional representation of the state of each cell that is corrected for batch effects
    and captures its gene expression profile and its environment.
-   Enables differential expression analysis across niches while accounting for wrong assignment of molecules to cells.
-   Scalable to very large datasets (>1 million cells).

The limitations of NicheVI include:

-   Effectively requires a GPU for fast inference.
-   Latent space is not interpretable, unlike that of a linear method.
-   Assumes single cells are observed and does not work with low resolution ST like Visium or Slide-Seq.

```{topic} Tutorials:

-   {doc}`/tutorials/notebooks/spatial/nicheVI_tutorial`
```

## Preliminaries

NicheVI takes as input spatially-resolved scRNA data. In addition to the gene expression matrix $\bm{X}$ with $N$ cells and $G$ genes,
it requires for each cell $n$:
- the spatial coordinates of the cell $y_n$
- the cell type assignment (possibly coarse) $c_n \in \{1, ..., T\}$
- the batch assignment $s_n$.


As preprocessing, we take the $K$ nearest neighbors of a cell to define its niche using the Euclidean distance in physical space.
We characterize the niche by its cell-type composition and gene expression. We denote by $\bm{\alpha_n}$ the $T$ dimensional vector of cell type
proportions among the $K$ nearest neighbors of the cell $n$. Its values are in the probability simplex.
The niche gene expression is defined as the average expression of each cell type present in the niche.
In practice, we leverage gene expression embeddings (PCA, scVI or similar) and characterize a cell type expression profile as the local average
embedding of cells of the same type. The average embeddings are stored in the matrix $\bm{\eta_n} \in \mathbb{R}^{T \times D}$, where $D$ is the embedding dimension.
## Descriptive model

We propose a latent variable model aiming to capture both gene expression heterogeneity and spatial variation resulting from the micro-environment.
We assume these two sources of variability are both captured by a $P$-dimensional latent variable ($P \ll G$):

```{math}
:nowrap: true
\begin{align}
    z_n \sim \mathbf{MixtureOfGaussians}(\mu_1, ..., \mu_M; \Sigma_1, ..,\Sigma_M; \pi_1, ...,\pi_M)
\end{align}
```

We model the observed expression profile $x_n$ as follows:

```{math}
:nowrap: true
\begin{align}
    x_n \mid z_n, s_n &\sim \mathbf{NegativeBinomial}\left( g_{\theta}(z_n, s_n)  \right)
\end{align}
```

The cell-type proportions of the cell's $K$ nearest neighbors are obtained as

```{math}
:nowrap: true
\begin{align}
    \alpha_n | z_n \sim \mathbf{Dirichlet}\left( g_{\omega}(z_n) \right),
\end{align}
```

Last, we assume that the neighboring cells' average expression profiles are obtained as

```{math}
:nowrap: true
\begin{equation}
\eta_{nt} | z_n, \alpha_n \sim \begin{cases}
\mathcal{N}\left(g_{\nu}^{t}(z_n)
      \right) , & \text{if $\alpha_{t} > 0$}\\
            0 & \text{otherwise}
    \end{cases}
\end{equation}
```
where $t=1,...,T$. $\theta$, $\omega$ and $\nu$ are neural network parameters.


The latent variables, along with their description are summarized in the following table:



## Inference

ResolVI uses variational inference, specifically auto-encoding variational Bayes
(see {doc}`/user_guide/background/variational_inference`) to learn both the model parameters
(the neural network parameters, dispersion parameters, etc.) and an approximate posterior distribution.

## Tasks

Here we provide an overview of some of the tasks that resolVI can perform. Please see {class}`scvi.external.RESOLVI`
for the full API reference.

### Dimensionality reduction

For dimensionality reduction, the mean of the approximate posterior $q_\phi(z_i \mid y_i, n_i)$ is returned by default.
This is achieved using the method:

```
>>> adata.obsm["X_resolvi"] = model.get_latent_representation()
```

Users may also return samples from this distribution, as opposed to the mean, by passing the argument `give_mean=False`.
The latent representation can be used to create a nearest neighbor graph with scanpy with:

```
>>> import scanpy as sc
>>> sc.pp.neighbors(adata, use_rep="X_resolvi")
>>> adata.obsp["distances"]
```

### Transfer learning

A resolVI model can be pre-trained on reference data and updated with query data using {meth}`~scvi.external.RESOLVI.load_query_data`, which then facilitates transfer of metadata like cell type annotations. $\beta_{N(n)n}$ is extended to the new cells and learned on these cells. The encoder by default does not see the batch covariate and $z_n$ can be predicted without performing query model training. See the {doc}`/user_guide/background/transfer_learning` guide for more information.

### Estimation of true expression levels

In {meth}`~scvi.external.RESOLVI.get_normalized_methylation` ResolVI returns the expected true expression value of $\rho_n$ under the approximate posterior. For one cell $n$, this can be written as:

```{math}
:nowrap: true

\begin{align}
   \mathbb{E}_{q_\phi(z_n \mid x_n)}\left[f_{\theta}\left(z_{n}, s_n \right) \right]
\end{align}
```

### Differential expression

Differential expression analysis is achieved with {meth}`~scvi.external.RESOLVI.differential_expression`.
ResolVI tests differences in expression levels $\rho_{n} = f_{\theta}\left(z_n, s_n\right)$.
We allow for importance based sampling using pyro's built-in function.

### Cell-type prediction

Prediction of cell-type labels is performed with {meth}`~scvi.external.RESOLVI.predict`.
A semisupervised model is necessary to perform this analysis as it leverages the cell-type classifier.
ResolVI performs for each cell $n$ $c_{n} = h_{nu}\left(z_n\right)$ and samples from $z_n$ to yield
the cell-type labels.

### Differential niche abundance

Differential niche abundance analysis is achieved with {meth}`~scvi.external.RESOLVI.differential_niche_abundance`.
A semisupervised model is necessary to perform this analysis as it leverages the cell-type classifier.
ResolVI tests differences in abundance of various cell-types in the neighborhood of a cell $n$
$c_{n} = h_{nu}\left(z_n\right)$. Cell-type prediction vectors are averaged weighted by the distance of a specific cell
and differential computation is performed.
