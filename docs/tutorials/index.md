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

:::{dropdown} Quick start tutorials

::::{grid} 1 2 3 3
:gutter: 2

:::{grid-item-card} [](notebooks/quick_start/api_overview)
:link: notebooks/quick_start/api_overview
:link-type: doc

Go through the typical steps of an scvi-tools workflow
:::

:::{grid-item-card} [](notebooks/quick_start/data_loading)
:link: notebooks/quick_start/data_loading
:link-type: doc

Load, preprocess, and register data for use with scvi-tools
:::

:::{grid-item-card} [](notebooks/quick_start/python_in_R)
:link: notebooks/quick_start/python_in_R
:link-type: doc

Perform basic Python operations in an R environment
:::

::::

:::

## scRNA-seq

Go through the workflows of models used for analysis of scRNA-seq datasets, including scVI, scANVI, AutoZI, AmortizedLDA, LinearSCVI, and contrastiveVI

### Pre-processing

### Unsupervised Models

:::{dropdown} Unsupervised model tutorials

:::

### Semi-supervised Models

:::{dropdown} semi-supervised model tutorials

:::

### Reference Mapping

:::{dropdown} Reference mapping tutorials

:::

## scATAC-seq

Go through the workflows of models used for analysis of scATAC-seq datasets, including PeakVI, scBasset, and PoissonVI

### Pre-processing

### PeakVI

:::{dropdown} PeakVI tutorials

:::

### PoissonVI

:::{dropdown} PoissonVI tutorials

:::

### scBasset

:::{dropdown} scBasset tutorials

:::

## scBS-seq (Single-cell bisulfite sequencing)

Use methylVI with single-cell bisulfite sequencing data

:::{dropdown} methylVI tutorials

:::

## Multimodal

Use models to analyze multimodal data, including totalVI for CITE-seq analysis and MultiVI for joint analysis of scRNA-seq and scATAC-seq

### Pre-processing

### totalVI

:::{dropdown} totalVI tutorials

:::

### MultiVI

:::{dropdown} MultiVI tutorials

:::

## Spatial Transcriptomics

Learn about models and methods to use with spatial transciptomics data, including DestVI, gimVI, Tangram, Cell2location, and Stereoscope

### Pre-processing

:::{dropdown} Spatial transciptomics tutorials

:::

## Model Hub

Learn about using and uploading pretrained scvi-tools models with scvi-hub and Hugging Face

:::{dropdown} Model hub tutorials

:::

## Common Modeling Use Cases

Learn common tasks in the scvi-tools workflow

:::{dropdown} Common use case tutorials

:::

## Development

Learn how to implement a novel statistical method for single-cell omics data in the scvi-tools environment

:::{dropdown} scvi-tools development tutorials

:::
