======
AUTOZI
======

**AUTOZI** [#ref1]_ (Python class :class:`scvi.model.AUTOZI`)
is a model for assessing gene-specific levels of zero-inflation in scRNA-seq data. 

Generative process
==================
AUTOZI is very similar to scVI but employs a spike-and-slab prior for the zero-inflation mixture assignment for each gene.
Whether the zero-inflation rate (:math:`\pi_{ng}` in the original scVI model) is sampled from a set of 
non-negligible values (the "slab" component) or the set of negligible values (the "spike" component) is defined by
:math:`m_g \sim Bernoulli(\delta_g)` where :math:`\delta_g \sim Beta(\alpha, \beta)`.
Thus, for the each gene :math:`g`, the zero-inflation rate is defined, 
:math:`\pi_{ng} = (1-m_g)\pi_{ng}^{slab} + m_g \pi_{ng}^{spike}`.

The full generative model is as follows:

.. math::
   :nowrap:

   \begin{align}
    z_n &\sim N(0,I)\\
    l_n &\sim LogNormal(l_u, l_\sigma^2)\\
    \delta_g &\sim Beta(\alpha^g,\beta^g)\\
    m_g &\sim Bernoulli(\delta_g)\\
    \pi _{ng} &=( 1-m_{g}) \delta _{\{0\}} +m_{g} \delta _{\{h^{g}( z_{n})\}}\\
    x_{ng}|z_n,l_n,m_g &\sim ZINB(l_nw_g(z_n), \theta_g, \pi_{ng})\\
    \end{align}

Where :math:`w^g` and :math:`h^g` are neural networks taking in :math:`z_n` and outputting 
the dropout rate and library size frequency respectively. The priors :math:`l_u` and 
:math:`l_{\sigma^2}` are the empircal mean and variance of the log library size per batch
respectively. The priors for :math:`\delta_g` are :math:`\alpha^g` and :math:`\beta^g` which 
by default are both set to 0.5 to enforce sparsity while maintaining symmetry. Finally,
:math:`\delta_{\{x\}}` denotes the Dirac distribution on :math:`x`.

Inference Procedure
===================

To learn the parameters, we employ variational inference with the following approximate posterior
distribution:

.. math::
   :nowrap:

    \begin{align*}
     \bar{q} &= \prod ^{G}_{g=1} q( \delta _{g})\prod ^{N}_{n=1} q( z_{n} |x_{n}) q( l_{n} |x_{n})
     \end{align*}

Tasks
=====
To classify whether a gene :math:`g` is or is not zero inflated, 
we call::

    >>> outputs = model.get_alpha_betas()
    >>> alpha_posterior = outputs['alpha_posterior']
    >>> beta_posterior = outputs['beta_posterior']

Then Bayesian decision theory suggests the posterior probability of of zero-inflation 
is :math:`q(\delta_g < 0.5)`.
    >>> from scipy.stats import beta
    >>> threshold = 0.5
    >>> zi_probs = beta.cdf(0.5, alpha_posterior, beta_posterior)