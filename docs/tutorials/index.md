# Tutorials

The easiest way to get familiar with scvi-tools is to follow along with our tutorials.
Many are also designed to work seamlessly in Google Colab, a free cloud computing platform.
Tutorials by default work with the latest installable version of scvi-tools. To view older tutorials,
change the documentation version using the tab at the bottom of the left sidebar.

:::{note}
For questions about using scvi-tools, or broader questions about modeling data, please use our [forum]. Checkout the [ecosystem] for additional models powered by scvi-tools.
:::

```{toctree}
:hidden:
:maxdepth: 2

index_quick_start
index_scrna
index_atac
index_scbs
index_multimodal
index_spatial
index_hub
index_use_cases
index_dev
```

[forum]: https://discourse.scverse.org/
[ecosystem]: https://scvi-tools.org/ecosystem


## Quick Start

Learn the typical scvi-tools workflow, how to handle data with scvi-tools, and how to use basic scvi-tools functionality within an R environment

::::{dropdown} Quick start tutorials

:::{card} [](notebooks/quick_start/api_overview)
:link: notebooks/quick_start/api_overview
:link-type: doc

Go through the typical steps of an scvi-tools workflow
:::

:::{card} [](notebooks/quick_start/data_loading)
:link: notebooks/quick_start/data_loading
:link-type: doc

Load, preprocess, and register data for use with scvi-tools
:::

:::{card} [](notebooks/quick_start/python_in_R)
:link: notebooks/quick_start/python_in_R
:link-type: doc

Perform basic Python operations in an R environment
:::

::::

## scRNA-seq

Go through the workflows of models used for analysis of scRNA-seq datasets, including scVI, scANVI, AutoZI, AmortizedLDA, LinearSCVI, and contrastiveVI

### Pre-processing

### Unsupervised Models

The following tutorials use generative models of scRNA-seq data for downstream analysis.
These models do not utilize known cell types.

::::{dropdown} Unsupervised model tutorials

:::{card}
:link: ./notebooks/scrna/harmonization
:link-type: doc

Perform integration of multiple scRNA-seq datasets both with and without cell type annotation (scVI and scANVI)
:::

:::{card}
:link: ./notebooks/scrna/scvi_in_R
:link-type: doc

Use basic scvi-tools functionality in R including integration of datasets
:::

:::{card}
:link: ./notebooks/scrna/scVI_DE_worm
:link-type: doc

Perform DE analysis on C. elegans data with scVI to quantify differences in gene expression between groups of cells
:::

:::{card}
:link: ./notebooks/scrna/MrVI_tutorial
:link-type: doc

Analyze multi-sample scRNA-seq data with MrVI
:::

:::{card}
:link: ./notebooks/scrna/sysVI
:link-type: doc

Integrate scRNA-seq datasets with substantial batch effects.
:::

:::{card}
:link: ./notebooks/scrna/contrastiveVI_tutorial
:link-type: doc

Use contrastiveVI to isolate perturbation-induced variation in Perturb-seq data
:::

:::{card}
:link: ./notebooks/scrna/amortized_lda
:link-type: doc

Run the amortized Latent Dirichlet Allocation model in scvi-tools to learn topics of an scRNA-seq dataset
:::

:::{card}
:link: ./notebooks/scrna/AutoZI_tutorial
:link-type: doc

Use the AutoZI model to enable gene-specific treatment of zero-inflation
:::

:::{card}
:link: ./notebooks/scrna/linear_decoder
:link-type: doc

Fit an LDVAE model to scRNA-seq data and interpret how genes are linked to latent variables of cells
:::

:::{card}
:link: ./notebooks/scrna/cellassign_tutorial
:link-type: doc

Use CellAssign to assign cell types using only knowledge of marker genes
:::

::::

### Semi-supervised Models

The following tutorials use scANVI, a semi-supervised model for single-cell transcriptomics data.
Unlike the models used in the previous section, scANVI can use existing cell type knowledge
for a subset of cells to infer the states of all cells.

::::{dropdown} semi-supervised model tutorials

:::{card} [](./notebooks/scrna/harmonization)
:link: ./notebooks/scrna/harmonization
:link-type: doc

Perform integration of multiple scRNA-seq datasets both with and without cell type annotation (scVI and scANVI)
:::

:::{card} [](./notebooks/scrna/seed_labeling)
:link: ./notebooks/scrna/seed_labeling
:link-type: doc

Create seed labels and transfer cell type annotations to an entire dataset
:::

:::{card} [](./notebooks/scrna/scanvi_fix)
:link: ./notebooks/scrna/scanvi_fix
:link-type: doc

Compare scANVI to other models following a bug fix in scvi-tools 1.1.0
:::

:::{card} [](./notebooks/scrna/tabula_muris)
:link: ./notebooks/scrna/tabula_muris
:link-type: doc

Perform de novo integration of a labeled reference dataset with an unlabeled query dataset (label transfer)
:::

::::

### Reference Mapping

The following tutorials cover the process of reference mapping with scArches. In the reference
mapping workflow, we have a model trained on reference data, and we map new query data to the
model's reference latent space in order to analyze the new data in the context of the reference
data. This is different from de novo integration, where the model is trained on multiple datasets.

::::{dropdown} Reference mapping tutorials

:::{card} [](./notebooks/scrna/scarches_scvi_tools)
:link: ./notebooks/scrna/scarches_scvi_tools
:link-type: doc

Map cells from a query dataset to the latent space of a reference dataset with the scArches method
:::

:::{card} [](./notebooks/scrna/query_hlca_knn)
:link: ./notebooks/scrna/query_hlca_knn
:link-type: doc

Use scANVI, scArches, and scvi-hub to query the Human Lung Cell Atlas
:::

::::

## scATAC-seq

Go through the workflows of models used for analysis of scATAC-seq datasets, including PeakVI, scBasset, and PoissonVI

### Pre-processing

### PeakVI

::::{dropdown} PeakVI tutorials

::::

### PoissonVI

::::{dropdown} PoissonVI tutorials

::::

### scBasset

::::{dropdown} scBasset tutorials

::::

## scBS-seq (Single-cell bisulfite sequencing)

Use methylVI with single-cell bisulfite sequencing data

::::{dropdown} methylVI tutorials

::::

## Multimodal

Use models to analyze multimodal data, including totalVI for CITE-seq analysis and MultiVI for joint analysis of scRNA-seq and scATAC-seq

### Pre-processing

### totalVI

::::{dropdown} totalVI tutorials

::::

### MultiVI

::::{dropdown} MultiVI tutorials

::::

## Spatial Transcriptomics

Learn about models and methods to use with spatial transciptomics data, including DestVI, gimVI, Tangram, Cell2location, and Stereoscope

### Pre-processing

::::{dropdown} Spatial transciptomics tutorials

::::

## Model Hub

Learn about using and uploading pretrained scvi-tools models with scvi-hub and Hugging Face

::::{dropdown} Model hub tutorials

::::

## Common Modeling Use Cases

Learn common tasks in the scvi-tools workflow

::::{dropdown} Common use case tutorials

::::

## Development

Learn how to implement a novel statistical method for single-cell omics data in the scvi-tools environment

::::{dropdown} scvi-tools development tutorials

::::
