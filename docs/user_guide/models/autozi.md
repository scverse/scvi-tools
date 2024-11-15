# AUTOZI

**AUTOZI** [^ref1] (Python class {class}`scvi.model.AUTOZI`)
is a model for assessing gene-specific levels of zero-inflation in scRNA-seq data.

```{topic} Tutorials:

-   {doc}`/tutorials/notebooks/scrna/AutoZI_tutorial`
```

## Generative process

AUTOZI is very similar to scVI but employs a spike-and-slab prior for the zero-inflation mixture assignment for each gene.
Whether the zero-inflation rate ($\pi_{ng}$ in the original scVI model) is sampled from a set of
non-negligible values (the "slab" component) or the set of negligible values (the "spike" component) is defined by
$m_g \sim Bernoulli(\delta_g)$ where $\delta_g \sim Beta(\alpha, \beta)$.
Thus, for each gene $g$, the zero-inflation rate is defined,
$\pi_{ng} = (1-m_g)\pi_{ng}^{slab} + m_g \pi_{ng}^{spike}$.

The full generative model is as follows:

```{math}
:nowrap: true

\begin{align}
 z_n &\sim N(0,I)\\
 l_n &\sim LogNormal(l_u, l_\sigma^2)\\
 \delta_g &\sim Beta(\alpha^g,\beta^g)\\
 m_g &\sim Bernoulli(\delta_g)\\
 \pi _{ng} &=( 1-m_{g}) \delta _{\{0\}} +m_{g} \delta _{\{h^{g}( z_{n})\}}\\
 x_{ng}|z_n,l_n,m_g &\sim ZINB(l_nw_g(z_n), \theta_g, \pi_{ng})\\
 \end{align}
```

Where $w^g$ and $h^g$ are neural networks taking in $z_n$ and outputting
the dropout rate and library size frequency respectively. The priors $l_u$ and
$l_{\sigma^2}$ are the empircal mean and variance of the log library size per batch
respectively. The priors for $\delta_g$ are $\alpha^g$ and $\beta^g$ which
by default are both set to 0.5 to enforce sparsity while maintaining symmetry. Finally,
$\delta_{\{x\}}$ denotes the Dirac distribution on $x$.

## Inference Procedure

To learn the parameters, we employ variational inference (see {doc}`/user_guide/background/variational_inference`) with the following approximate posterior
distribution:

```{math}
:nowrap: true

 \begin{align*}
  \bar{q} &= \prod ^{G}_{g=1} q( \delta _{g})\prod ^{N}_{n=1} q( z_{n} |x_{n}) q( l_{n} |x_{n})
  \end{align*}
```

## Tasks

To classify whether a gene $g$ is or is not zero inflated,
we call:

```
>>> outputs = model.get_alpha_betas()
>>> alpha_posterior = outputs['alpha_posterior']
>>> beta_posterior = outputs['beta_posterior']
```

Then Bayesian decision theory suggests the posterior probability of of zero-inflation
is $q(\delta_g < 0.5)$.

```
>>> from scipy.stats import beta
>>> threshold = 0.5
>>> zi_probs = beta.cdf(0.5, alpha_posterior, beta_posterior)
```

[^ref1]: Oscar Clivio, Romain Lopez, Jeffrey Regier, Adam Gayoso, Michael I. Jordan, Nir Yosef (2019), _Detecting zero-inflated genes in single-cell transcriptomics data_, [Machine Learning in Computational Biology (MLCB)](https://www.biorxiv.org/content/biorxiv/early/2019/10/10/794875.full.pdf).
