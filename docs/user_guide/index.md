# User guide

scvi-tools is composed of models that can perform one or many analysis tasks. In the user guide, we provide an overview of each model with emphasis on the math behind the model, how it connects to the code, and how the code connects to analysis.

:::{figure} /\_static/tasks.png
:align: center
:alt: Overview of tasks
:class: img-fluid
:::

## scRNA-seq analysis

```{eval-rst}
.. list-table::
   :widths: 15 15 100
   :header-rows: 1

   * - Model
     - Reference
     - Tasks
   * - :doc:`/user_guide/models/scvi`
     - :cite:p:`Lopez18`
     - Dimensionality reduction, removal of unwanted variation, integration across replicates, donors, and technologies, differential expression, imputation, normalization of other cell- and sample-level confounding factors
   * - :doc:`/user_guide/models/scanvi`
     - :cite:p:`Xu21`
     - scVI tasks with cell type transfer from reference, seed labeling
   * - :doc:`/user_guide/models/linearscvi`
     - :cite:p:`Svensson20`
     - scVI tasks with linear decoder
   * - :doc:`/user_guide/models/autozi`
     - :cite:p:`Clivio19`
     -  for assessing gene-specific levels of zero-inflation in scRNA-seq data
   * - :doc:`/user_guide/models/cellassign`
     - :cite:p:`Zhang19`
     - Marker-based automated annotation
   * - :doc:`/user_guide/models/solo`
     - :cite:p:`Bernstein20`
     - Doublet detection
   * - :doc:`/user_guide/models/scar`
     - :cite:p:`Sheng22`
     - Ambient RNA removal

```

## ATAC-seq analysis

```{eval-rst}
.. list-table::
   :widths: 15 15 100
   :header-rows: 1

   * - Model
     - Reference
     - Tasks
   * - :doc:`/user_guide/models/peakvi`
     - :cite:p:`Ashuach22`
     - Dimensionality reduction, removal of unwanted variation, integration across replicates, donors, and technologies, differential expression, imputation, normalization of other cell- and sample-level confounding factors
   * - :doc:`/user_guide/models/scbasset`
     - :cite:p:`Yuan2022`
     - Dimensionality reduction, removal of unwanted variation, integration across replicates, donors, and technologies, imputation
```

## Multimodal analysis

### CITE-seq

```{eval-rst}
.. list-table::
   :widths: 15 15 100
   :header-rows: 1

   * - Model
     - Reference
     - Tasks
   * - :doc:`/user_guide/models/totalvi`
     - :cite:p:`GayosoSteier21`
     - Dimensionality reduction, removal of unwanted variation, integration across replicates, donors, and technologies, differential expression, protein imputation, imputation, normalization of other cell- and sample-level confounding factors
```

### Multiome

```{eval-rst}
.. list-table::
   :widths: 15 15 100
   :header-rows: 1

   * - Model
     - Reference
     - Tasks
   * - :doc:`/user_guide/models/multivi`
     - :cite:p:`AshuachGabitto21`
     - Integration of paired/unpaired multiome data, missing modality imputation, normalization of other cell- and sample-level confounding factors

```

## Spatial transcriptomics analysis

```{eval-rst}
.. list-table::
   :widths: 15 15 100
   :header-rows: 1

   * - Model
     - Reference
     - Tasks
   * - :doc:`/user_guide/models/destvi`
     - :cite:p:`Lopez21`
     - Multi-resolution deconvolution, cell-type-specific gene expression imputation, comparative analysis
   * - :doc:`/user_guide/models/stereoscope`
     - :cite:p:`Andersson20`
     - Deconvolution
   * - :doc:`/user_guide/models/gimvi`
     - :cite:p:`Lopez19`
     - Imputation of missing spatial genes
   * - :doc:`/user_guide/models/tangram`
     - :cite:p:`Biancalani21`
     - Deconvolution, single cell spatial mapping
```

## General purpose analysis

```{eval-rst}
.. list-table::
   :widths: 15 15 100
   :header-rows: 1

   * - Model
     - Reference
     - Tasks
   * - :doc:`/user_guide/models/amortizedlda`
     - :cite:p:`Blei03`
     - Topic modeling

```

## Background

-   {doc}`/user_guide/background/variational_inference`
-   {doc}`/user_guide/background/differential_expression`
-   {doc}`/user_guide/background/counterfactual_prediction`
-   {doc}`/user_guide/background/transfer_learning`
-   {doc}`/user_guide/background/codebase_overview`
