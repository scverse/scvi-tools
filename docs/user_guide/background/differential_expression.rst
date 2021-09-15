==============================
Differential Expression
==============================

Under construction.

Problem statement
==========================================

Differential expression analyses aim to quantify and detect expression differences of some quantity between conditions, e.g., cell types.
In single-cell experiments, such quantity can correspond to transcripts, protein expression, or chromatin accessibility.
A central notion when comparing expression levels of two cell states 
is the log fold-change

.. math::
   :nowrap:

   \begin{align}
      \beta_g := \log h_{g}^B - \log h_{g}^A,
   \end{align}

where 
:math:`\log h_{g}^A, \log h_{g}^B`
respectively denote the mean expression levels in subpopulations :math:`A`
and
:math:`B`.



Motivations to use scVI-tools for differential expression 
======================================================================

In the particular case of single-cell RNA-seq data, existing differential expression models often model that the mean expression level 
:math:`\log h_{g}^C`.
as linear functions of the cell-state and batch assignments.
These models face two notable limitations to detect differences in expression between cell-states in large-scale scRNA-seq datasets.
First, such linear assumptions may not capture complex batch effects existing in such datasets accurately.
When comparing two given states :math:`A`
and
:math:`B` in a large dataset, these models may also struggle to leverage data of a related state present in the data.

Deep generative models may not suffer from these problems.
Most ``scvi-tools`` models use complex nonlinear mappings to capture batch effects on expression.
Using amortization, they can leverage large amounts of data
to better capture shared correlations between features.
Consequently, deep generative models have appealing properties for differential expression in large-scale data.

This guide has two objectives.
First, it aims to provide insight as to how scVI-tools' differential expression module works for transcript expression (``scVI``), surface protein expression (``TOTALVI``), or chromatin accessibility (``PeakVI``).
More precisely, we explain how it can:

    + approximate population-specific normalized expression levels

    + detect biologically relevant features

    + provide easy-to-interpret predictions

More importantly, this guide explains the function of the hyperparameters of the ``differential_expression`` method.


.. list-table::
   :widths: 20 50 15 15 15
   :header-rows: 1

   * - Parameter
     - Description
     - Approximating expression levels
     - Detecting relevant features
     - Providing easy-to-interpret predictions
   * - ``idx1``, ``idx2``
     - Mask or queries for the compared populations :math:`A` and :math:`B`.
     - yes
     - 
     - 
   * - ``mode``
     - Characterizes the null hypothesis.
     - 
     - yes
     - 
   * - ``delta``
     - composite hypothesis characteristics (when ``mode="change"``).
     - 
     - yes
     - 
   * - ``fdr_target``
     - desired FDR significance level
     - 
     - 
     - yes
   * - ``importance_sampling``
     - Precises if expression levels are estimated using importance sampling
     - yes
     - 
     - 

Notations and model assumptions
======================================================================
While considering different modalities, scVI, TOTALVI, and PeakVI share similar properties, allowing us to perform differential expression of transcripts, surface proteins, or chromatin accessibility, similarly.
We first introduce some notations that will be useful in the remainder of this guide.

[COMPLETE]

 

Approximating population-specific normalized expression levels
====================================================================================

A first step to characterize differences in expression consists in estimating state-specific expression levels.
Most ``scVI-tools`` models do not explicitly model discrete cell types for several reasons. 
The most obvious one is that states often are unknown at the beginning of the analysis.
In some cases, states may have an intricate structure that would be difficult to model.
Instead, the class of models we consider here assumes that a latent variable :math:`z` characterizes cells' biological identity.
A key component of our differential expression module is to aggregate the information carried by individual cells to estimate population-wide expression levels.
The strategy to do so is as follows.
First, we estimate latent representations of the two compared states :math:`A` and :math:`B`, using aggregate variational posteriors.
In particular, we will represent state :math:`C` latent representation with the mixture

.. math::
   :nowrap:

   \begin{align}
      \hat P^C(
        Z
      ) = 
      \frac
      {1}
      {
        \mathcal{N}_C
      }
      \sum_{n \in \mathcal{N}_C}
      p_\theta(z \mid x_n),
   \end{align}

``idx1`` and``idx2`` specify which observations to use to approximate these quantities.

Detecting biologically relevant features
========================================================


Providing easy-to-interpret predictions
========================================================


Understanding the differential expression output
========================================================
