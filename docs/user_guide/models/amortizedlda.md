# Amortized LDA

**LDA** [^ref1] (Latent Dirichlet Allocation) posits a generative model where
a set of latent topics generates collections of elements. In the case of single-cell RNA sequencing, we can think
of these topics as gene modules and each cell as a collection of UMI counts. Other features that can be ascribed to these
topics include surface proteins and accessible chromatin regions, all of which have discrete count values.
This implementation (Python class {class}`~scvi.model.AmortizedLDA`) of LDA amortizes the
cost of performing variational inference for each cell by training a common encoder. Note: this is not an exact implementation
of the model described in the original LDA paper.

The advantages of amortized LDA are:

-   Can learn underlying topics without a reference.
-   Scalable to very large datasets (>1 million cells).

The limitations of amortized LDA include:

-   Optimal selection of the number of topics is unclear.
-   Amortization gap in optimizing variational parameters.

```{topic} Tutorials:

-   {doc}`/tutorials/notebooks/scrna/amortized_lda`
```

## Preliminaries

Amortized LDA takes as input a cell-by-feature matrix $X$ with $C$ cells and $F$ features.
Because the LDA model assumes the input is ordered, we refer to this format as the bag-of-words (BoW) representation
of the feature counts.
Additionally, the number of topics to model must be manually set by the user prior to fitting the model.

## Generative process

Amortized LDA posits that the $N$ observed feature counts for cell $c$ are treated as ordered. For all $n \in [N]$ feature counts
for cell $c \in [C]$, the observed feature counts $x_{cn}$ are produced according to the following generative process:

```{math}
:nowrap: true

\begin{align}
 \beta_k &\sim \mathrm{Dir}(\eta)  &\forall k \in [K]\\
 \theta_c &\sim \mathrm{Dir}(\alpha) &\\
 x_{cn} &\sim \mathrm{Cat}(\theta_c \beta) &\forall n \in [N] \\
\end{align}
```

where $\eta$ denotes the prior on the Dirichlet distribution for the topic feature distribution $\beta$,
and $\alpha$ denotes the prior on the Dirichlet distribution for the cell topic distribution $\theta_c$.
In order to compute reparametrization gradients stably, we approximate the Dirichlet distribution with a logistic-Normal
distribution, followed by a softmax operation. Specifically, we use the Laplace approximation
which has a diagonal covariance matrix [^ref2]:

```{math}
:nowrap: true

\begin{align}
 \mu_k &= \log\alpha_k - \frac{1}{K}\sum_i \log\alpha_i \\
 \Sigma_{kk} &= \frac{1}{\alpha_k} \left(1 - \frac{2}{K}\right) + \frac{1}{K^2} \sum_i \frac{1}{\alpha_k}
\end{align}
```

for Dirichlet parameter $\alpha \in \mathbb{R}^K$ where $K$ denotes the number of topics.

The latent variables, along with their description are summarized in the following table:

```{eval-rst}
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
     - Parameter for the Dirichlet prior on the topic feature distribution, :math:`\beta_k`. Approximated by a logistic-Normal distribution.
     - ``topic_feature_prior``
   * - :math:`\theta_c \in \Delta^{K-1}`
     - Cell topic distribution for a given cell :math:`c`.
     - ``cell_topic_dist``
   * - :math:`\beta_k \in \Delta^{F-1}`
     - Topic feature distribution for a given topic :math:`k`.
     - ``topic_feature_dist``
```

## Inference

Amortized LDA uses variational inference and specifically auto-encoding variational bayes (see {doc}`/user_guide/background/variational_inference`)
to learn both the model parameters (the neural network params, topic feature distributions, etc.) and an approximate posterior distribution.
Like {class}`scvi.model.SCVI`, the underlying class used as the encoder for Amortized LDA is {class}`~scvi.nn.Encoder`.

## Tasks

### Topic-based dimensionality reduction

Users can retrieve the estimated topic proportions in each cell with the following code:

```
>>> topic_prop = model.get_latent_representation()
>>> adata.obsm["X_LDA"] = topic_prop
```

Due to the logistic-Normal distribution not having an analytic solution to the mean, we compute
a Monte Carlo estimate of the expectation. The number of samples used for the estimate can be configured
with the argument `n_samples`.

Additionally, once can estimate topic proportions on held-out data by passing in an AnnData object
with the same format as the dataset used to train the model:

```
>>> test_topic_prop = model.get_latent_representation(test_adata)
```

If the learned topics generalize well to other datasets, this can serve as a dimensionality reduction method
to the learned topic latent space.

### Feautre module discovery

Once the model has been fitted, one can retrieve the estimated feature-by-topic distribution:

```
>>> feature_by_topic = model.get_feature_by_topic()
```

Like the `get_latent_representation()` method, this returns a Monte Carlo estimate of the logistic-Normal expectation.
Similarly, we can configure the number of samples with `n_samples`.

[^ref1]:
    David M. Blei, Andrew Y. Ng, Michael I. Jordan (2003),
    _Latent Dirichlet Allocation_,
    [Journal of Machine Learning Research](https://www.jmlr.org/papers/volume3/blei03a/blei03a.pdf).

[^ref2]:
    Akash Srivastava, Charles Sutton (2017),
    _Autoencoding Variational Inference for Topic Models_,
    [International Conference on Learning Representations](https://arxiv.org/pdf/1703.01488.pdf).
