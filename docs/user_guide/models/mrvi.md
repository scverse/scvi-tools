# MrVI

**MrVI** [^ref1] (Multi-resolution Variational Inference; Python class
{class}`~scvi.external.MRVI`) is a deep generative model designed for the analysis of large-scale
single-cell transcriptomics data with multi-sample, multi-batch experimental designs.

MrVI conducts both **exploratory analyses** (locally dividing samples into groups based on molecular properties)
and **comparative analyses** (comparing pre-defined groups of samples in terms of differential expression and differential abundance) at single-cell resolution.
It can capture nonlinear and cell-type specific variation of sample-level covariates on gene expression.

```{topic} Tutorials:

-    {doc}`/tutorials/notebooks/scrna/MrVI_tutorial`
```

## Preliminaries

MrVI takes as input a scRNA-seq gene expression matrix $X$ with $N$ cells and $G$ genes.
Additionally, it requires specification, for each cell $n$:
- a sample-level target covariate $s_n$, that typically corresponds to the sample ID,
	which defines which sample entities will be compared in exploratory and comparative analyses.
- nuisance covariates $b_n$ (e.g. sequencing run, processing day).

Optionally, MrVI can also take as input
	- Cell-type labels for guided integration across samples, via a mixture of Gaussians prior where each mixture component serves as the mode of a cell type.
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

:::{figure} figures/mrvi_graphical_model.svg
:align: center
:alt: MrVI graphical model
:class: img-fluid

MrVI graphical model. Shaded nodes represent observed data, unshaded nodes represent latent variables.
:::

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
MrVI also enables supervised comparative analysis to detect cell-type specific DE and DA between sample groups.

#### Differential expression
At a high level, the DE procedure regresses, within each cell $n$, counterfactual cell states $z^{(s)}_n$ on sample-level covariates $c_s$ of interest for analysis as
$z^{(s)}_n = \beta_n c_s + \beta_0 + \epsilon_n$.
For instance, if $c_s$ is a binary covariate, then $\beta_n$ will capture the shift (in $z$-space) induced by samples for which $c_s = 1$ compared to samples for which $c_s = 0$.
This procedure, repeated for all cells, allows several downstream analyses.
First, comparing the norm of $\beta_n$ (using $\chi^2$ statistics) across cells can identify cell-states that vary the most for a given covariate, or conversely, identify sample covariates that strongly associate with specific cell states.
Second, by decoding the linear approximation of $z^{(s)}_n$ for different covariate vectors that we would like to compare, we can compute associated log fold-changes to identify DE genes at the cell level.

#### Differential abundance
To compare two sets of samples, MrVI computes the log-ratio between the aggregated posteriors of the two groups, $A_1 \subset [[1, S]]$ and $A_2 \subset [[1, S]]$, where $S$ is the total number of samples.
In particular, the aggregated posterior for any sample $s$ is defined as
$q_s := \frac{1}{|s|} \sum_{n, s_n=s} q^{u}_{n}$,
where $q_n$ is the posterior of cell $n$ in $u$-space.
This aggregated posterior $q_s$ characterizes the distribution of all cells in sample $s$.
To characterize the distribution of cells in a group of samples $A$, we can consider the mixture of aggregated posteriors $q_s$ for all $s \in A$, corresponding to
$q_A := \frac{1}{|A|} \sum_{s \in A} q_s$.
In particular, cell states $u$ that are abundant in a sample group $A$ will have a high probability mass in $q_A$, while rare states will have low probability mass.
More generally, we can consider the log-ratio of aggregated posteriors between two groups of samples $A_1$ and $A_2$ as a measure of differential abundance:
$r = \log \frac{q_{A_1}}{q_{A_2}}$.
We can evaluate these log-ratios for all cell states $u$ to identify DA cell-state regions.
In particular, large positive (resp. negative) values of $r$ indicate that cell states are more abundant in $A_1$ (resp. $A_2$).

[^ref1]:
     Pierre Boyeau, Justin Hong, Adam Gayoso, Martin Kim, Jose L McFaline-Figueroa, Michael Jordan, Elham Azizi, Can Ergen, Nir Yosef (2024),
    _Deep generative modeling of sample-level heterogeneity in single-cell genomics_,
    [bioRxiv](https://doi.org/10.1101/2022.10.04.510898).
