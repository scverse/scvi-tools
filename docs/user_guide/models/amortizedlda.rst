======
Amortized LDA
======

**LDA** [#ref1]_ (Latent Dirichlet Allocation; Python class :class:`~scvi.model.AmortizedLDA`) posits a generative model where
a set of latent topics generates collections of elements. In the case of single-cell RNA sequencing, we can think
of these topics as gene modules and each cell as a collection of UMI counts. This implementation of LDA amortizes the
cost of performing variational inference for each cell by training a common encoder.

The advantages of amortized LDA are:

    + Can learn underlying topics without a reference.

    + Scalable to very large datasets (>1 million cells).

The limitations of amortized LDA include:

    + Optimal selection of the number of topics is unclear.

    + Amortization gap in optimizing variational parameters.

.. topic:: Tutorials:

 - :doc:`/tutorials/notebooks/amortized_lda`


Preliminaries
==============
Amortized LDA takes as input a scRNA-seq gene expression matrix :math:`X` with :math:`N` cells and :math:`G` genes.
Because the LDA model assumes the input is ordered, we refer to this format as the bag-of-words (BoW) representation
of the gene expression UMI counts.
Additionally, the number of topics to model must be manually set by the user prior to fitting the model.


Generative process
==================

Amortized LDA posits that the :math:`N` observed UMI counts for cell :math:`c` are treated as ordered. For the :math:`n`th UMI count,
the gene expressed, :math:`x_{cn}` is produced according to the following generative process:

.. math::
   :nowrap:

   \begin{align}
    \text{For each topic $k$:} \\
      &\beta_k \sim \mathrm{Dir}(\eta) \\
    \theta_c &\sim \mathrm{Dir}(\alpha) \\
    \text{For each UMI count $n$:} \\
      &x_{cn} \sim \mathrm{Cat}(\theta_c \beta_{z_{cn}})
   \end{align}

where :math:`\eta` denotes the prior on the Dirichlet distribution for the topic gene distribution :math:`\beta`,
and :math:`\alpha` denotes the prior on the Dirichlet distribution for the cell topic distribution :math:`\theta_c`.
In order to compute reparametrization gradients stably, we approximate the Dirichlet distribution with a logistic-Normal
distribution, followed by a softmax operation. Specifically, we use the Laplace approximation
which has a diagonal covariance matrix [#ref2]_:

.. math::
   :nowrap:

   \begin{align}
    \mu_k &= \mathrm{log}\alpha_k - \frac{1}{K}\sum_i \mathrm{log} \alpha_i \\
    \Sigma_{kk} &= \frac{1}{\alpha_k} \left(1 - \frac{2}{K}\right) + \frac{1}{K^2} \sum_i \frac{1}{\alpha_k}
   \end{align}

for Dirichlet parameter :math:`\alpha \in \mathbb{R}^K` where :math:`K` denotes the number of topics.

The latent variables, along with their description are summarized in the following table:

.. list-table::
   :widths: 20 90 15
   :header-rows: 1

   * - Latent variable
     - Description
     - Code variable (if different)
   * - :math:`\alpha \in (0, \infty)^K`
     - Parameter for the Dirichlet prior on the cell topic distribution, :math:`\theta_c`. Approximated by a logistic-Normal distribution.
     - ``cell_topic_prior``
   * - :math:`\eta \in (0, \infty)^K`
     - Parameter for the Dirichlet prior on the topic gene distribution, :math:`\beta_k`. Approximated by a logistic-Normal distribution.
     - ``topic_gene_prior``
   * - :math:`\theta_c \in \Delta^{K-1}`
     - Cell topic distribution for a given cell :math:`c`.
     - ``cell_topic_dist``
   * - :math:`\beta_k \in \Delta^{G-1}`
     - Topic gene distribution for a given topic :math:`k`.
     - ``topic_gene_dist``

Inference
=========

Amortized LDA uses variational inference and specifically auto-encoding variational bayes (see :doc:`/user_guide/background/variational_inference`)
to learn both the model parameters (the neural network params, topic gene distributions, etc.) and an approximate posterior distribution.
Like :class:`scvi.model.SCVI`, the underlying class used as the encoder for Amortized LDA is :class:`~scvi.nn.Encoder`.

Tasks
=====

Topic-based dimensionality reduction
------------------------------------

Users can retrieve the estimated topic proportions in each cell with the following code:

    >>> topic_prop = model.get_latent_representation()
    >>> adata.obsm["X_LDA"] = topic_prop

Due to the logistic-Normal distribution not having an analytic solution to the mean, we compute
a Monte Carlo estimate of the expectation. The number of samples used for the estimate can be configured
with the argument ``n_samples``.

Additionally, once can estimate topic proportions on held-out data by passing in an AnnData object
with the same format as the dataset used to train the model:

    >>> test_topic_prop = model.get_latent_representation(test_adata)

If the learned topics generalize well to other datasets, this can serve as a dimensionality reduction method
to the learned topic latent space.

Gene module discovery
---------------------

Once the model has been fitted, one can retrieve the estimated gene-by-topic distribution:

    >>> gene_by_topic = model.get_gene_by_topic()

Like the ``get_latent_representation()`` method, this returns a Monte Carlo estimate of the logistic-Normal expectation.
Similarly, we can configure the number of samples with ``n_samples``.

.. topic:: References:

   .. [#ref1] David M. Blei, Andrew Y. Ng, Michael I. Jordan (2003),
      *Latent Dirichlet Allocation*,
      `Journal of Machine Learning Research <https://www.jmlr.org/papers/volume3/blei03a/blei03a.pdf>`__.
   .. [#ref2] Akash Srivastava, Charles Sutton (2017),
      *Autoencoding Variational Inference for Topic Models*,
      `International Conference on Learning Representations <https://arxiv.org/pdf/1703.01488.pdf>`__.
