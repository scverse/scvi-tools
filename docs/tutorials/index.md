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

:::{card} []()
:link: ./notebooks/scrna/harmonization
:link-type: doc

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

The following tutorials use peakVI, a generative model of scATAC-seq data.

::::{dropdown} PeakVI tutorials

:::{card} [](./notebooks/atac/PeakVI)
:link: ./notebooks/atac/PeakVI
:link-type: doc

Go through the PeakVI workflow to analyze a scATAC-seq dataset
:::

:::{card} [](./notebooks/atac/peakvi_in_R)
:link: ./notebooks/atac/peakvi_in_R
:link-type: doc

Use scvi-tools functionality in R to analyze scATAC-seq data
:::

::::

### PoissonVI

In the following tutorial, we use poissonVI, which models scATAC-seq fragment counts, unlike
peakVI which only models binarized scATAC-seq data.

::::{dropdown} PoissonVI tutorials

:::{card} [](./notebooks/atac/PoissonVI)
:link: ./notebooks/atac/PoissonVI
:link-type: doc

Go through the PoissonVI workflow to analyze scATAC-seq data using quantitative fragment counts
:::

::::

### scBasset

In the following tutorials, we use scBasset which provides a sequence-based method for
representation learning of scATAC-seq data.

::::{dropdown} scBasset tutorials

:::{card} [](./notebooks/atac/scbasset)
:link: ./notebooks/atac/scbasset
:link-type: doc

Go through the scBasset workflow to analyze a scATAC-seq dataset
:::

:::{card} [](./notebooks/atac/scbasset_batch)
:link: ./notebooks/atac/scbasset_batch
:link-type: doc

Use scBasset to integrate data across several samples
:::

::::

## scBS-seq (Single-cell bisulfite sequencing)

Use methylVI with single-cell bisulfite sequencing data

::::{dropdown} methylVI tutorials

:::{card} [](notebooks/scbs/MethylVI_batch)
:link: notebooks/scbs/MethylVI_batch
:link-type: doc

Correct batch effects in across different scBS-seq experiments with methylVI
:::

::::

## Multimodal

Use models to analyze multimodal data, including totalVI for CITE-seq analysis and MultiVI for joint analysis of scRNA-seq and scATAC-seq

### Pre-processing

### totalVI

In the following tutorials we use totalVI, a generative model of CITE-seq RNA and protein data
that can be used for downstream analysis tasks such as dimensionality reduction and differential
expression. TotalVI can also be used to integrate CITE-seq and scRNA-seq data.

::::{dropdown} totalVI tutorials

:::{card} [](./notebooks/multimodal/totalVI)
:link: ./notebooks/multimodal/totalVI
:link-type: doc

Go through the totalVI workflow to analyze CITE-seq datasets
:::

:::{card} [](./notebooks/multimodal/totalVI_reference_mapping)
:link: ./notebooks/multimodal/totalVI_reference_mapping
:link-type: doc

Use totalVI to train a reference model and map CITE-seq query data
:::

:::{card} [](./notebooks/multimodal/cite_scrna_integration_w_totalVI)
:link: ./notebooks/multimodal/cite_scrna_integration_w_totalVI
:link-type: doc

Use totalVI to integrate CITE-seq and scRNA-seq datasets
:::

:::{card} [](./notebooks/multimodal/totalvi_in_R)
:link: ./notebooks/multimodal/totalvi_in_R
:link-type: doc

Use scvi-tools functionality in R to analyze CITE-seq data
:::

::::

### MultiVI

In the following tutorial, we use MultiVI, a multimodal generative model which can integrate
multiome, scRNA-seq and scATAC-seq data. MultiVI can be used for downstream tasks such as
dimensionality reduction and imputation of missing modalities.

::::{dropdown} MultiVI tutorials

:::{card} [](./notebooks/multimodal/MultiVI_tutorial)
:link: ./notebooks/multimodal/MultiVI_tutorial
:link-type: doc

Go through the MultiVI workflow to perform joint analysis of paired and unpaired multi omic data
:::

::::

## Spatial Transcriptomics

Learn about models and methods to use with spatial transciptomics data, including DestVI, gimVI, Tangram, Cell2location, and Stereoscope

### Pre-processing

::::{dropdown} Spatial transciptomics tutorials

:::{card}
In the following tutorials we use various scvi-tools models for spatial transcriptomics
analysis. See the [user guide](./user_guide/index.md#Spatial transcriptomics analysis) for more information about these models.
:::

:::{card} [](./notebooks/spatial/resolVI_tutorial)
:link: ./notebooks/spatial/resolVI_tutorial
:link-type: doc

description
:::

:::{card} [](./notebooks/spatial/DestVI_tutorial)
:link: ./notebooks/spatial/DestVI_tutorial
:link-type: doc

Perform multi-resolution analysis on spatial transcriptomics data with DestVI
:::

:::{card} [](./notebooks/spatial/DestVI_in_R)
:link: ./notebooks/spatial/DestVI_in_R
:link-type: doc

Use scvi-tools functionality in R to analyze spatial transcriptomics datasets
:::

:::{card} [](./notebooks/spatial/gimvi_tutorial)
:link: ./notebooks/spatial/gimvi_tutorial
:link-type: doc

Use gimVI to impute missing genes in spatial data
:::

:::{card} [](./notebooks/spatial/tangram_scvi_tools)
:link: ./notebooks/spatial/tangram_scvi_tools
:link-type: doc

Use Tangram to map spatial transcriptomics data
:::

:::{card} [](./notebooks/spatial/stereoscope_heart_LV_tutorial)
:link: ./notebooks/spatial/stereoscope_heart_LV_tutorial
:link-type: doc

Go through the Stereoscope workflow to map single-cell data
:::

:::{card} [](./notebooks/spatial/cell2location_lymph_node_spatial_tutorial)
:link: ./notebooks/spatial/cell2location_lymph_node_spatial_tutorial
:link-type: doc

Spatially map lymph node cell types using Cell2location
:::

::::

## Model Hub

Learn about using and uploading pretrained scvi-tools models with scvi-hub and Hugging Face

::::{dropdown} Model hub tutorials

:::{card} [](./notebooks/hub/scvi_hub_intro_and_download)
:link: ./notebooks/hub/scvi_hub_intro_and_download
:link-type: doc

Learn how to use Hugging Face and scvi-hub to download pretrained scvi-tools models
:::

:::{card} [](./notebooks/hub/cellxgene_census_model)
:link: ./notebooks/hub/cellxgene_census_model
:link-type: doc

Perform analysis of a CELLxGENE dataset using a pretrained model from scVI-hub
:::

:::{card} [](./notebooks/hub/scvi_hub_upload_and_large_files)
:link: ./notebooks/hub/scvi_hub_upload_and_large_files
:link-type: doc

Learn how to upload pretrained scvi-tools models to Hugging Face
:::

::::

## Common Modeling Use Cases

Learn common tasks in the scvi-tools workflow

::::{dropdown} Common use case tutorials

:::{card} [](./notebooks/use_cases/autotune_scvi)
:link: ./notebooks/use_cases/autotune_scvi
:link-type: doc

Automatically find a good set of hyperparameters using autotune.
:::

:::{card} [](./notebooks/hub/minification)
:link: ./notebooks/hub/minification
:link-type: doc

Minify a dataset by replacing count data with the model’s estimated parameters of the latent posterior distribution
:::

:::{card} [](./notebooks/use_cases/interpretability)
:link: ./notebooks/use_cases/interpretability
:link-type: doc

Interpret a model to improve its performance
:::

:::{card} [](./notebooks/use_cases/external_datasets)
:link: ./notebooks/use_cases/external_datasets
:link-type: doc

Learn how to use scvi-tools models for external datasets
:::

::::

## Development

Learn how to implement a novel statistical method for single-cell omics data in the scvi-tools environment

::::{dropdown} scvi-tools development tutorials

:::{card} [](./notebooks/dev/data_tutorial)
:link: ./notebooks/dev/data_tutorial
:link-type: doc

Learn about how data is handled in scvi-tools
:::

:::{card} [](./notebooks/dev/module_user_guide)
:link: ./notebooks/dev/module_user_guide
:link-type: doc

Implement a novel statistical method for single-cell omics data as a module
:::

:::{card} [](./notebooks/dev/model_user_guide)
:link: ./notebooks/dev/model_user_guide
:link-type: doc

Implement an scvi-tools model class to provide a convenient interface for the lower-level module objects
:::

::::
