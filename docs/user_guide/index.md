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
   :widths: 15 100 25
   :header-rows: 1

   * - Model
     - Tasks
     - Reference
   * - :doc:`/user_guide/models/scvi`
     - Dimensionality reduction, removal of unwanted variation, integration across replicates, donors, and technologies, differential expression, imputation, normalization of other cell- and sample-level confounding factors
     - :cite:p:`Lopez18`
   * - :doc:`/user_guide/models/scanvi`
     - scVI tasks with cell type transfer from reference, seed labeling
     - :cite:p:`Xu21`
   * - :doc:`/user_guide/models/linearscvi`
     - scVI tasks with linear decoder
     - :cite:p:`Svensson20`
   * - :doc:`/user_guide/models/autozi`
     -  for assessing gene-specific levels of zero-inflation in scRNA-seq data
     - :cite:p:`Clivio19`
   * - :doc:`/user_guide/models/cellassign`
     - Marker-based automated annotation
     - :cite:p:`Zhang19`
   * - :doc:`/user_guide/models/solo`
     - Doublet detection
     - :cite:p:`Bernstein20`
   * - :doc:`/user_guide/models/scar`
     - Ambient RNA removal
     - :cite:p:`Sheng22`
   * - :doc:`/user_guide/models/contrastivevi`
     - scVI tasks with contrastive analysis
     - :cite:p:`Weinberger23`
   * - :doc:`/user_guide/models/mrvi`
     - Characterization of sample-level heterogeneity
     - :cite:p:`Boyeau24`
   * - :doc:`/user_guide/models/sysvi`
     - Integrating single-cell RNA-seq datasets with substantial batch effects
     - :cite:p:`Hrovatin23`
   * - :doc:`/user_guide/models/decipher`
     - Joint representation and visualization of derailed cell states with Decipher
     - :cite:p:`Nazaret24`
   * - :doc:`/user_guide/models/velovi`
     - Deep generative modeling of transcriptional dynamics for RNA velocity analysis in single cells
     - :cite:p:`GayosoWeiler23`
```

## ATAC-seq analysis

```{eval-rst}
.. list-table::
   :widths: 15 100 25
   :header-rows: 1

   * - Model
     - Tasks
     - Reference
   * - :doc:`/user_guide/models/peakvi`
     - Dimensionality reduction, removal of unwanted variation, integration across replicates, donors, and technologies, differential expression, imputation, normalization of other cell- and sample-level confounding factors
     - :cite:p:`Ashuach22`
   * - :doc:`/user_guide/models/scbasset`
     - Dimensionality reduction, removal of unwanted variation, integration across replicates, donors, and technologies, imputation
     - :cite:p:`Yuan2022`
   * - :doc:`/user_guide/models/poissonvi`
     - Dimensionality reduction, removal of unwanted variation, integration across replicates, donors, and technologies, differential expression, imputation, normalization of other cell- and sample-level confounding factors
     - :cite:p:`Martens2023`
```

## BS-seq analysis

```{eval-rst}
.. list-table::
   :widths: 15 100 25
   :header-rows: 1

   * - Model
     - Tasks
     - Reference
   * - :doc:`/user_guide/models/methylvi`
     - Dimensionality reduction, removal of unwanted variation, integration across replicates, donors, and technologies, differential methylation, imputation, normalization of other cell- and sample-level confounding factors
     - :cite:p:`Weinberger2023a`
   * - :doc:`/user_guide/models/methylanvi`
     - MethylVI tasks along with cell type label transfer from reference, seed labeling
     - :cite:p:`Weinberger2023a`
```

## Multimodal analysis

### CITE-seq

```{eval-rst}
.. list-table::
   :widths: 15 100 25
   :header-rows: 1

   * - Model
     - Tasks
     - Reference
   * - :doc:`/user_guide/models/totalvi`
     - Dimensionality reduction, removal of unwanted variation, integration across replicates, donors, and technologies, differential expression, protein imputation, imputation, normalization of other cell- and sample-level confounding factors
     - :cite:p:`GayosoSteier21`
```

### Multiome

```{eval-rst}
.. list-table::
   :widths: 15 100 25
   :header-rows: 1

   * - Model
     - Tasks
     - Reference
   * - :doc:`/user_guide/models/multivi`
     - Integration of paired/unpaired multiome data, missing modality imputation, normalization of other cell- and sample-level confounding factors
     - :cite:p:`AshuachGabitto21`

```

## Spatial transcriptomics analysis

```{eval-rst}
.. list-table::
   :widths: 15 100 25
   :header-rows: 1

   * - Model
     - Tasks
     - Reference
   * - :doc:`/user_guide/models/destvi`
     - Multi-resolution deconvolution, cell-type-specific gene expression imputation, comparative analysis
     - :cite:p:`Lopez22`
   * - :doc:`/user_guide/models/stereoscope`
     - Deconvolution
     - :cite:p:`Andersson20`
   * - :doc:`/user_guide/models/gimvi`
     - Imputation of missing spatial genes
     - :cite:p:`Lopez19`
   * - :doc:`/user_guide/models/tangram`
     - Deconvolution, single cell spatial mapping
     - :cite:p:`Biancalani21`
   * - :doc:`/user_guide/models/resolvi`
     - Generative model of single-cell resolved spatial transcriptomics
     - :cite:p:`Ergen25`
```

## General purpose analysis

```{eval-rst}
.. list-table::
   :widths: 15 100 25
   :header-rows: 1

   * - Model
     - Tasks
     - Reference
   * - :doc:`/user_guide/models/amortizedlda`
     - Topic modeling
     - :cite:p:`Blei03`

```

## Background

- {doc}`/user_guide/background/variational_inference`
- {doc}`/user_guide/background/differential_expression`
- {doc}`/user_guide/background/counterfactual_prediction`
- {doc}`/user_guide/background/transfer_learning`
- {doc}`/user_guide/background/codebase_overview`

## Common Use Cases

- {doc}`/user_guide/use_case/saving_and_loading_models`
- {doc}`/user_guide/use_case/downstream_analysis_tasks`
- {doc}`/user_guide/use_case/using_callbacks`
- {doc}`/user_guide/use_case/hyper_parameters_tuning`
- {doc}`/user_guide/use_case/multi_gpu_training`
- {doc}`/user_guide/use_case/custom_dataloaders`
- {doc}`/user_guide/use_case/scvi_criticism`

## Glossary

::::{tab-set}

:::{tab-item} Model

A Model class inherits {class}`~scvi.model.base.BaseModelClass` and is the user-facing object for interacting with a module.
The model has a `train` method that learns the parameters of the module, and also contains methods
for users to retrieve information from the module, like the latent representation of cells in a VAE.
Conventionally, the post-inference model methods should not store data into the AnnData object, but
instead return "standard" Python objects, like numpy arrays or pandas dataframes.
:::

:::{tab-item} Module

A module is the lower-level object that defines a generative model and inference scheme. A module will
either inherit {class}`~scvi.module.base.BaseModuleClass` or {class}`~scvi.module.base.PyroBaseModuleClass`.
Consequently, a module can either be implemented with PyTorch alone, or Pyro. In the PyTorch only case, the
generative process and inference scheme are implemented respectively in the `generative` and `inference` methods,
while the `loss` method computes the loss, e.g, ELBO in the case of variational inference.
:::
::::

::::{tab-set}

:::{tab-item} TrainingPlan

The training plan is a PyTorch Lightning Module that is initialized with a scvi-tools module object.
It configures the optimizers, defines the training step and validation step, and computes metrics to be
recorded during training. The training step and validation step are functions that take data, run it through
the model and return the loss, which will then be used to optimize the model parameters in the Trainer.
Overall, custom training plans can be used to develop complex inference schemes on top of modules.
:::

:::{tab-item} Trainer

The {class}`~scvi.train.Trainer` is a lightweight wrapper of the PyTorch Lightning Trainer. It takes as input
the training plan, a training data loader, and a validation dataloader. It performs the actual training loop, in
which parameters are optimized, as well as the validation loop to monitor metrics. It automatically handles moving
data to the correct device (CPU/GPU).
:::
::::
