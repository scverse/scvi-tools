# MrVI

**MrVI** [^ref1] (Multi-resolution Variational Inference; Python class
{class}`~scvi.external.MRVI`) is a deep generative model designed for the analysis of large-scale
single-cell transcriptomics data with multi-sample, multi-batch experimental designs.

The advantages of MrVI are:

- Facilitates both exploratory analysis to divide samples into groups based on
    molecular properties and comparative analysis to compare pre-defined groups of samples.
- Provides single-cell-resolution differential expression and differential abundance estimates.
- Captures nonlinear, cell-type specific effects of sample-level covariates on gene expression.
- Scales to large datasets with millions of cells and hundreds of samples.

The limitations of MrVI include:

- Requires specification of sample-level target (i.e., of biological interest) and nuisance (i.e., introduces undesired variation) covariates.
- Differential expression procedure does not provide single-cell gene-level p-values.
MrVI conducts both **exploratory analyses** (locally dividing samples into groups based on molecular properties)
and **comparative analyses** (comparing pre-defined groups of samples in terms of differential expression and differential abundance) at single-cell resolution.
It can capture nonlinear and cell-type specific variation of sample-level covariates on gene expression.
```{topic} Tutorials:

-    {doc}`tutorials/notebooks/scrna/MrVI_tutorial`
```

## Preliminaries

MrVI takes as input a scRNA-seq gene expression matrix $X$ with $N$ cells and $G$ genes.
Additionally, it requires specification, for each cell $n$:
- a sample-level target covariate $s_n$, that typically corresponds to the sample ID,
	which defines which sample entities will be compared in explorary and comparative analyses.
- nuisance covariates $b_n$ (e.g. sequencing run, processing day).

Optionally, MrVI can also take as input
	- Cell-type labels for improved integration across samples, via a guided mixture of Gaussians prior.
	- Additional sample-level covariates of interest $c_s$ for each sample $s$ (e.g.
	  disease status, age, treatment) for comparative analysis.

## Generative process

MrVI posits a two-level hierarchical model (Figure 1):

1. A cell-level latent variable $u_n$ capturing cell state in a batch-corrected manner:
    $u_n \sim \mathrm{MixtureOfGaussians}(\mu_1, ..., \mu_K, \Sigma_1, ..., \Sigma_K, \pi_1, ..., \pi_K)$
2. A cell-level latent variable $z_n$ capturing both cell state and effects of sample $s_n$:
    $z_n | u_n \sim \mathcal{N}(u_n, I_L)$
3. The normalized gene expression levels $h_n$ are generated from $z_n$ as:
    $h_n = \mathrm{softmax}(A_{zh} \times [z_n + g_\theta(z_n, b_n)] + \gamma_{zh})$
4. Finally the gene expression counts are generated as:
    $x_{ng} | h_{ng} \sim \mathrm{NegativeBinomial}(l_n h_{ng}, r_{ng})$

Here $l_n$ is the library size of cell $n$, $r_{ng}$ is the gene-specific inverse dispersion,
$A_{zh}$ is a linear matrix of dimension $G \times L$, $\gamma_{zh}$ is a bias vector of dimension
$G$, and $\theta$ are neural network parameters.
$u_n$ captures broad cell states invariant to sample and batch,
while $z_n$ augments $u_n$ with sample-specific effects while correcting for nuisance covariate effects.
Gene expression is obtained from $z_n$ using multi-head attention mechanisms to
    flexibly model batch and sample effects.

The latent variables, along with their description are summarized in the following table:

```{eval-rst}
.. list-table::
   :widths: 20 90 15
   :header-rows: 1

   * - Latent variable
     - Description
     - Code variable (if different)
   * - :math:`u_n \in \mathbb{R}^L`
  - "sample-unaware" representation of a cell, invariant to both sample and nuisance covariates.
     - ``u``
   * - :math:`z_n \in \mathbb{R}^L`
  - "sample-aware" representation of a cell, invariant to nuisance covariates.
     - ``z``
   * - :math:`h_n \in \mathbb{R}^G`
     - Cell-specific normalized gene expression.
     - ``h``
   * - :math:`l_n \in \mathbb{R}^+`
     - Cell size factor.
     - ``library``
   * - :math:`r_{ng} \in \mathbb{R}^+`
     - Gene and cell-specific inverse dispersion.
     - ``px_r``
   * - :math:`\mu_1, ..., \mu_K \in \mathbb{R}^L`
     - Mixture of Gaussians means for prior on $u_n$.
     - ``u_prior_means``
   * - :math:`\Sigma_1, ..., \Sigma_K \in \mathbb{R}^{L \times L}`
     - Mixture of Gaussians covariance matrices for prior on $u_n$.
     - ``u_prior_scales``
   * - :math:`\pi_1, ..., \pi_K \in \mathbb{R}^+`
     - Mixture of Gaussians weights for prior on $u_n$.
     - ``u_prior_logits``
```

## Inference

MrVI uses variational inference to approximate the posterior of $u_n$ and $z_n$. The variational
distributions are:

$q_{\phi}(u_n | x_n) := \mathcal{N}(\mu_{\phi}(x_n), \sigma^2_{\phi}(x_n)I)$

$z_n := u_n + f_{\phi}(u_n, s_n)$

Here $\mu_{\phi}, \sigma^2_{\phi}$ are encoder neural networks and $f_{\phi}$ is a deterministic
mapping based on multi-head attention between $u_n$ and a learned embedding for sample $s_n$.

## Tasks

### Exploratory analysis

MrVI enables unsupervised local sample stratification via the construction of cell-specific
sample-sample distance matrices, for every cell $n$:

1. For each cell state $u_n$, compute counterfactual cell states $z^{(s)}_n$ for all possible samples $s$.
2. Compute cell-specific sample-sample distance matrices $D^{(n)}$ based on the Euclidean distance between all pairs of $z^{(s)}_n$.
3. Cluster cells based on their $D^{(n)}$ to find cell populations with distinct sample stratifications.
4. Average $D^{(n)}$ within each cell cluster and hierarchically cluster samples
This automatically reveals distinct sample stratifications that are specific to particular cell
subsets.

### Comparative analysis

MrVI also enables supervised comparative analysis to detect cell-type specific DE and DA between
sample groups:

The differential expression procedure at a high level is:

1. For each cell $n$, regress $z^{(s)}_n$ on sample-level covariates $c_s$ to obtain cell-level
    coefficients $\beta_n$,
2. Identify significant covariates based on $\chi^2$ statistic of $\beta_n$,
3. Predict and decode $z^{(s')}_n$ for given covariate vectors and compute log fold-changes to get
    cell-level DE genes.

The differential abundance procedure at a high level is:

1. Compute aggregated posteriors for each sample $q_s(u)$,
2. For each sample group, average the densities of the constituent samples to get a
    covariate-group-specific aggregated posterior,
3. Compare sample groups by computing the log probability ratios between the
    covariate-group-specific aggregated posteriors.

[^ref1]:
    TODO: Add reference here
