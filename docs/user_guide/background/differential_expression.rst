==============================
Differential Expression
==============================

Under construction.

Problem statement
==========================================

Differential expression consists in quantifying and detecting gene expression differences between conditions, e.g. cell-types.
To do so, a central notion when comparing expression levels of two cell states 
:math:`A`
and
:math:`B` is the log fold-change

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



Motivations to use scVI for differential expression 
======================================================================

Existing differential expression models often model that the mean expression level 
:math:`\log h_{g}^C`.
as linear functions of the cell-state and batch assignments.
These models face two notable limitations to detect differences in expression between cell-states in large-scale scRNA-seq datasets.
First, such linear assumptions may not properly capture complex batch effects existing in such datasets.
When comparing two given states :math:`A`
and
:math:`B` in a large dataset, these models may also struggle to leverage data of a related state 
:math:`C`.
present in the data.

Deep generative models often do not suffer from these problems.
In scVI, batch effects on gene expression are modelled nonlinearly by a neural network.
Consequently, scVI has appealing properties to detect differentially expressed genes in large scale data.

This guide has two objectives.
First, it aims to provide insight as of how scVI's differential expression module works, and more precisely how it was designed to:

    + approximate population-specific normalized expression levels

    + detect biologically relevant genes

    + provide easy-to-interpret genes

More importantly, this guide explains the function of the hyperparameters of the ``differential_expression`` method.


.. list-table::
   :widths: 20 50 15 15 15
   :header-rows: 1

   * - Parameter
     - Description
     - Approximating expression levels
     - Detecting relevant genes
     - Providing easy-to-interpret genes
   * - ``idx1``, ``idx2``
     - N/A
     - yes
     - 
     - 
   * - ``mode``
     - N/A
     - 
     - yes
     - 
   * - ``delta``
     - N/A
     - 
     - yes
     - 
   * - ``fdr_target``
     - N/A
     - 
     - 
     - yes
   * - ``importance_sampling``
     - N/A
     - yes
     - 
     - 
   * - ``pseudocounts``
     - N/A
     - 
     - yes
     - 

Approximating population-specific normalized expression levels
====================================================================================

scVI does not explicitly model discrete cell-types. 
Instead, it assumes that a latent variable $z$ contains the biological

Detecting biologically relevant genes
========================================================


Providing easy-to-interpret genes
========================================================


Understanding scVI's differential expression output
========================================================
