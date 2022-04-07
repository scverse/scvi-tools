# User guide

scvi-tools is composed of models that can perform one or many analysis tasks. In the user guide, we provide an overview of each model with emphasis on the math behind the model, how it connects to the code, and how the code connects to analysis.

:::{figure} /_static/tasks.png
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
     - [Lopez18]_
     - Dimensionality reduction, removal of unwanted variation, integration across replicates, donors, and technologies, differential expression, imputation, normalization of other cell- and sample-level confounding factors
   * - :doc:`/user_guide/models/scanvi`
     - [Xu21]_
     - scVI tasks with cell type transfer from reference, seed labeling
   * - :doc:`/user_guide/models/linearscvi`
     - [Svensson20]_
     - scVI tasks with linear decoder
   * - :doc:`/user_guide/models/autozi`
     - [Clivio19]_
     -  for assessing gene-specific levels of zero-inflation in scRNA-seq data
   * - :doc:`/user_guide/models/cellassign`
     - [Zhang19]_
     - Marker-based automated annotation
   * - :doc:`/user_guide/models/solo`
     - [Bernstein19]_
     - Doublet detection

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
     - [Ashuach22]_
     - Dimensionality reduction, removal of unwanted variation, integration across replicates, donors, and technologies, differential expression, imputation, normalization of other cell- and sample-level confounding factors
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
     - [GayosoSteier21]_
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
     - [AshuachGabitto21]_
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
     - [Lopez21]_
     - Multi-resolution deconvolution, cell-type-specific gene expression imputation, comparative analysis
   * - :doc:`/user_guide/models/stereoscope`
     - [Andersson20]_
     - Deconvolution
   * - :doc:`/user_guide/models/gimvi`
     - [Lopez19]_
     - Imputation of missing spatial genes
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
     - [Blei03]_
     - Topic modeling

```

## Background

- {doc}`/user_guide/background/variational_inference`
- {doc}`/user_guide/background/differential_expression`
- {doc}`/user_guide/background/counterfactual_prediction`
- {doc}`/user_guide/background/transfer_learning`
- {doc}`/user_guide/background/codebase_overview`
