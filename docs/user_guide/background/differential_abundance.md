# Differential Abundance

:::{note}
This page is under construction.
:::

## Problem Statement

Differential abundance analyses aim to detect cell states that are disproportionally abundant in a given group of samples.
In the differential abundance procedure used in scvi-tools, we estimate the posterior density for each sample in the sample-unaware latent space (e.g. the u latent space in MrVI). These estimated densities can then be evaluated at particular cells, providing useful information about the abundance of specific cell states in various samples. In particular, for two disjoint sets of samples $A$ and $B$, we quantify the relative overabundance of cells from $A$ compared with $B$ at any cell state $u$ by computing the log of the aggregated posterior densities of the two groups:

```{math}
:nowrap: true

\begin{align}
   r_{AB}(u) := \log{\frac{q_A(u)}{q_B(u)}}
\end{align}

```

## Motivation

More importantly, this guide explains the function of the parameters of the `differential_abundance` method.

```{eval-rst}
.. list-table::
   :widths: 20 50 15 15 15
   :header-rows: 1

    * - Parameter
      - Description
      - Other
      - Other
    * - sample_key
      - Key for the sample covariate.
      -
      -
    * - num_cells_posterior
      - Maximum possible number of cells used to compute aggregated posterior for each sample. This should be used to avoid running out of memory if very large samples are present, as the aggregated posteriors are costly and can easily cause out of memory errors.
      -
      -
    * - dof
      - Degrees of freedom for the Student's t-distribution components for aggregated posterior. If ``None``, components are Normal. Using a Student's t-distribution instead of a Normal distribution to estimate posteriors can improve later analysis by creating a smoother distribution for samples containing few cells (in this case our estimated posterior comes from aggregating across fewer cells, so the estimate might be more innacurate, particularly between cells).
      -
      -

```

## Notations and model assumptions

While different `scvi-tools` models may consider different modalities (gene expression, protein expression, multimodal, etc.), they share similar properties, namely some low-dimensional representation of each cell. In particular, we consider a deep generative model where a latent variable with prior $z_n$ represents cell $n$'s identity. In turn, a neural network $f^h_\theta$ maps this low-dimensional representation to normalized, expression levels.

## Quantifying the probability that a given cell belongs a given sample

The first step to identifying abundant cell states in a given group of samples consists of, for each cell $n$ in a given sample $s$, using the model's variational approximation to the posterior distribution for $n$ over the $z$ space, $q(z|x_n)$. Aggregating these posteriors over all cells in $s$ results in the following aggregated posterior distribution:

```{math}
:nowrap: true

\begin{align}
   q_s(z) := \frac{1}{n_s}\sum_{n:s_n = s}q(z|x_n)
\end{align}

```
where $n_s$ is the number of cells in $s$. By evaluating this density function at the $z$ space representation of cell $n$, $q_s(z_n)$, we obtain the probability that cell $n$ belongs to sample $s$. In `scvi-tools`, the differential abundance method results in a $log$ probabilities matrix containing the $logs$ of these probabilities between each cell and sample.

## Aggregating posteriors to identify relatively overabundant cell states in a given group of samples

Next, we can quantify the density of any set of samples $A$ in the $z$ space as $q_A(z) := \frac{1}{|A|} \sum_{s \in A}q_s(z)$. Now, by evaluating this density function at the $z$ space representation of a given cell $n$, we obtain the probability that cell $n$ belongs to the group of samples, $A$. In practice this is useful because we can use covariates such as sex, or patient condition, to group the samples. By assigning groups this way, we identify cell states that are relatively overabundant in certain covariate groups, such as diseased individuals or females.

## Sources
https://docs.scvi-tools.org/en/1.4.1/user_guide/models/mrvi.html
https://docs.scvi-tools.org/en/1.4.1/user_guide/models/cytovi.html
