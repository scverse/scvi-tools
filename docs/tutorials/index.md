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

## scRNA-seq

Go through the workflows of models used for analysis of scRNA-seq datasets, including scVI, scANVI, AutoZI, AmortizedLDA, LinearSCVI, and contrastiveVI

### Pre-processing

### Unsupervised Models

### Semi-supervised Models

### Reference Mapping

## scATAC-seq

Go through the workflows of models used for analysis of scATAC-seq datasets, including PeakVI, scBasset, and PoissonVI

### Pre-processing

### PeakVI

### PoissonVI

### scBasset

## scBS-seq (Single-cell bisulfite sequencing)

Use methylVI with single-cell bisulfite sequencing data

## Multimodal

Use models to analyze multimodal data, including totalVI for CITE-seq analysis and MultiVI for joint analysis of scRNA-seq and scATAC-seq

### Pre-processing

### totalVI

### MultiVI

## Spatial Transcriptomics

Learn about models and methods to use with spatial transciptomics data, including DestVI, gimVI, Tangram, Cell2location, and Stereoscope

### Pre-processing

## Model Hub

Learn about using and uploading pretrained scvi-tools models with scvi-hub and Hugging Face

## Common Modeling Use Cases

Learn common tasks in the scvi-tools workflow

## Development

Learn how to implement a novel statistical method for single-cell omics data in the scvi-tools environment
