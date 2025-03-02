# Tutorials

The easiest way to get familiar with scvi-tools is to follow along with our tutorials.
Many are also designed to work seamlessly in Google Colab, a free cloud computing platform.
Tutorials by default work with the latest installable version of scvi-tools. To view older tutorials,
change the documentation version using the tab at the bottom of the left sidebar.

:::{note}
For questions about using scvi-tools, or broader questions about modeling data, please use our [forum]. Checkout the [ecosystem] for additional models powered by scvi-tools.
:::

```{toctree}
:maxdepth: 2
:hidden:

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

:::{dropdown} Quick start


Learn the typical scvi-tools workflow, how to handle data with scvi-tools, and how to use basic scvi-tools functionality within an R environment

:::{list-table}
:widths: auto
:align: left

*   - [](./notebooks/quick_start/api_overview)
    - Go through the typical steps of an scvi-tools workflow
*   - [](./notebooks/quick_start/data_loading)
    - Load, preprocess, and register data for use with scvi-tools
*   - [](./notebooks/quick_start/python_in_R)
    - Perform basic Python operations in an R environment

:::

:::

:::{dropdown} scRNA-seq

Go through the workflows of models used for analysis of scRNA-seq datasets, including scVI, scANVI, AutoZI, AmortizedLDA, LinearSCVI, and contrastiveVI

:::{card} Unsupervised Models
The following tutorials use generative models of scRNA-seq data for downstream analysis.
These models do not utilize known cell types.
:::

:::{list-table}
:widths: auto
:align: left

*   - Preprocessing
    - description
*   - [](./notebooks/scrna/harmonization) (link to first header)
    - Perform integration of multiple scRNA-seq datasets both with and without cell type annotation (scVI and scANVI)
*   - [](./notebooks/scrna/scvi_in_R)
    - Use basic scvi-tools functionality in R including integration of datasets
*   - [](./notebooks/scrna/scVI_DE_worm)
    - Perform DE analysis on C. elegans data with scVI to quantify differences in gene expression between groups of cells
*   - [](./notebooks/scrna/MrVI_tutorial)
    - Analyze multi-sample scRNA-seq data with MrVI
*   - [](./notebooks/scrna/MrVI_tutorial)
    - Preprocess data and train the MrVI model
*   - [](./notebooks/scrna/MrVI_tutorial)
    - Visualize MrVI’s u latent space
*   - [](./notebooks/scrna/MrVI_tutorial)
    - Perform differential expression and differential abundance analysis
*   - [](./notebooks/scrna/contrastiveVI_tutorial)
    - Use contrastiveVI to isolate perturbation-induced variation in Perturb-seq data
*   - [](./notebooks/scrna/amortized_lda)
    - Run the amortized Latent Dirichlet Allocation model in scvi-tools to learn topics of an scRNA-seq dataset
*   - [](./notebooks/scrna/AutoZI_tutorial)
    - Use the AutoZI model to enable gene-specific treatment of zero-inflation
*   - [](./notebooks/scrna/linear_decoder)
    - Fit an LDVAE model to scRNA-seq data and interpret how genes are linked to latent variables of cells
*   - [](./notebooks/scrna/cellassign_tutorial)
    - Use CellAssign to assign cell types using only knowledge of marker genes

:::

:::{card} scANVI
The following tutorials use scANVI, a semi-supervised model for single-cell transciptomics data.
Unlike the models used in the previous section, scANVI can use exisiting cell type knowledge
for a subset of cells to infer the states of all cells.
:::

:::{list-table}

*   - Preprocessing
    - description
*   - [](./notebooks/scrna/harmonization) (link to second header)
    - Perform integration of multiple scRNA-seq datasets both with and without cell type annotation (scVI and scANVI)
*   - [](./notebooks/scrna/seed_labeling)
    - Create seed labels and transfer cell type annotations to an entire dataset
*   - [](./notebooks/scrna/scanvi_fix)
    - Compare scANVI to other models following a bug fix in scvi-tools 1.1.0
*   - [](./notebooks/scrna/tabula_muris)
    - Perform de novo integration of a labeled reference dataset with an unlabeled query dataset (label transfer)

:::

:::{card} Reference Mapping
The following tutorials cover the process of reference mapping with scArches. In the reference
mapping workflow, we have a model trained on reference data, and we map new query data to the
model's reference latent space in order to analyze the new data in the context of the reference
data. This is different from de novo integration, where the model is trained on multiple datasets.
:::

:::{list-table}

*   - [](./notebooks/scrna/scarches_scvi_tools)
    - Map cells from a query dataset to the latent space of a reference dataset with the scArches method
*   - [](./notebooks/scrna/query_hlca_knn)
    - Use scANVI, scArches and scvi-hub to query the Human Lung Cell Atlas

:::
:::

:::{dropdown} scATAC-seq



Go through the workflows of models used for analysis of scATAC-seq datasets, including PeakVI, scBasset, and PoissonVI



:::{card}
The following tutorials use peakVI, a generative model of scATAC-seq data.
:::

:::{list-table}
:widths: auto
:align: left

*   - Preprocessing
    - description
*   - [](./notebooks/atac/PeakVI)
    - Go through the PeakVI workflow to analyze a scATAC-seq dataset
*   - [](./notebooks/atac/peakvi_in_R)
    - Use scvi-tools functionality in R to analyze scATAC-seq data

:::

:::{card}
In the following tutorial, we use poissonVI, which models scATAC-seq fragment counts, unlike
peakVI which only models binarized scATAC-seq data.
:::

:::{list-table}
:widths: auto
:align: left

*   - [](./notebooks/atac/PoissonVI)
    - Go through the PoissonVI workflow to analyze scATAC-seq data using quantitative fragment counts

:::

:::{card}
In the following tutorials, we use scBasset which provides a sequence-based method for
representation learning of scATAC-seq data.
:::

:::{list-table}
:widths: auto
:align: left

*   - [](./notebooks/atac/scbasset)
    - Go through the scBasset workflow to analyze a scATAC-seq dataset
*   - [](./notebooks/atac/scbasset_batch)
    - Use scBasset to integrate data across several samples

:::
:::

:::{dropdown} scBS-seq (Single-cell bisulfite sequencing)



Use methylVI with single-cell bisulfite sequencing data

:::{list-table}
:widths: auto
:align: left

*   - [](notebooks/scbs/MethylVI_batch)
    - Correct batch effects in across different scBS-seq experiments with methylVI
*   - [](notebooks/scbs/MethylVI_batch)
    - Correct batch effects in across different scBS-seq experiments with methylVI

:::
:::

:::{dropdown} Multimodal



Use models to analyze multimodal data, including totalVI for CITE-seq analysis and MultiVI for joint analysis of scRNA-seq and scATAC-seq



:::{card}
In the following tutorials we use totalVI, a generative model of CITE-seq RNA and protein data
that can be used for downstream analysis tasks such as dimensionalit reduction and differential
expression. TotalVI can also be used to integrate CITE-seq and scRNA-seq data.
:::

:::{list-table}
:widths: auto
:align: left

*   - Preprocessing
    - description
*   - [](./notebooks/multimodal/totalVI)
    - Go through the totalVI workflow to analyze CITE-seq datasets
*   - [](./notebooks/multimodal/totalVI_reference_mapping)
    - Use totalVI to train a reference model and map CITE-seq query data
*   - [](./notebooks/multimodal/cite_scrna_integration_w_totalVI)
    - Use totalVI to integrate CITE-seq and scRNA-seq datasets
*   - [](./notebooks/multimodal/totalvi_in_R)
    - Use scvi-tools functionality in R to analyze CITE-seq data

:::

:::{card}
In the following tutorial, we use MultiVI, a multimodal generative model which can integrate
multiome, scRNA-seq and scATAC-seq data. MultiVI can be used for downstream tasks such as
dimensionality reduction and imputation of missing modalities.
:::

:::{list-table}
:widths: auto
:align: left

*   - Preprocessing
    - description
*   - [](./notebooks/multimodal/MultiVI_tutorial)
    - Go through the MultiVI workflow to perform joint analysis of paired and unpaired multi omic data

:::
:::

:::{dropdown} Spatial transcriptomics


Learn about models and methods to use with spatial transciptomics data, including DestVI, gimVI, Tangram, Cell2location, and Stereoscope



:::{card}
In the following tutorials we use various scvi-tools models for spatial transcriptomics
analysis. See the [user guide](./user_guide/index.md#Spatial transcriptomics analysis) for more information about these models.
:::

:::{list-table}
:widths: auto
:align: left

*   - Preprocessing
    - description
*   - [](./notebooks/spatial/resolVI_tutorial)
    - description
*   - [](./notebooks/spatial/DestVI_tutorial)
    - Perform multi-resolution analysis on spatial transcriptomics data with DestVI
*   - [](./notebooks/spatial/DestVI_in_R)
    - Use scvi-tools functionality in R to analyze spatial transcriptomics datasets
*   - [](./notebooks/spatial/gimvi_tutorial)
    - Use gimVI to impute missing genes in spatial data
*   - [](./notebooks/spatial/tangram_scvi_tools)
    - Use Tangram to map spatial transcriptomics data
*   - [](./notebooks/spatial/stereoscope_heart_LV_tutorial)
    - Go through the Stereoscope workflow to map single-cell data
*   - [](./notebooks/spatial/cell2location_lymph_node_spatial_tutorial)
    - Spatially map lymph node cell types using Cell2location

:::
:::

:::{dropdown} Model hub



Learn about using and uploading pretrained scvi-tools models with scvi-hub and Hugging Face

:::{list-table}
:widths: auto
:align: left


*   - [](./notebooks/hub/minification)
    - Minify a dataset by replacing count data with the model’s estimated parameters of the latent posterior distribution
*   - [](./notebooks/hub/scvi_hub_intro_and_download)
    - Learn how to use Hugging Face and scvi-hub to download pretrained scvi-tools models
*   - [](./notebooks/hub/cellxgene_census_model)
    - Perform analysis of a CELLxGENE dataset using a pretrained model from scVI-hub
*   - [](./notebooks/hub/scvi_hub_upload_and_large_files)
    - Learn how to upload pretrained scvi-tools models to Hugging Face



:::
:::

:::{dropdown} Hyperparameter tuning



Learn about hyperparameter tuning of scvi-tools models

:::{list-table}
:widths: auto
:align: left

*   - [](./notebooks/tuning/autotune_scvi)
    - Use scvi’s autotune module to automatically find a good set of model hyperparameters

:::
:::

:::{dropdown} Development

Learn how to implement a novel statistical method for single-cell omics data in the scvi-tools environment

:::{list-table}
:widths: auto
:align: left

*   - [](./notebooks/dev/data_tutorial)
    - Learn about how data is handled in scvi-tools
*   - [](./notebooks/dev/module_user_guide)
    - Implement a novel statistical method for single-cell omics data as a module
*   - [](./notebooks/dev/model_user_guide)
    - Implement an scvi-tools model class to provide a convenient interface for the lower-level module objects

:::
:::

[forum]: https://discourse.scverse.org/
[ecosystem]: https://scvi-tools.org/ecosystem
