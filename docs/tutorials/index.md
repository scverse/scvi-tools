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
index_tuning
index_dev
```

:::{list-table} Tutorial Tags
:widths: auto

*   - <span class="tag" group-tag="Analysis"></span> Analysis
    - <span class="tag" group-tag="Integration"></span> Integration
    - <span class="tag" group-tag="Transfer learning"></span> Transfer learning
    - <span class="tag" group-tag="R"></span> R
    - <span class="tag" group-tag="Differential comparison"></span> Differential comparison
*   - <span class="tag" group-tag="Modality imputation"></span> Modality imputation
    - <span class="tag" group-tag="Deconvolution"></span> Deconvolution
    - <span class="tag" group-tag="Dimensionality reduction"></span> Dimensionality reduction
    - <span class="tag" group-tag="Removal of unwanted variance"></span> Removal of unwanted variance
    -
:::

:::{dropdown} Quick start

Learn the typical scvi-tools workflow, how to handle data with scvi-tools, and how to use basic scvi-tools functionality within an R environment

:::{list-table}
:widths: auto
:align: left

*   - <span class="tag" group-tag="Analysis"></span>
    <span class="tag" group-tag="Differential comparison"></span>
    <span class="tag" group-tag="Dimensionality reduction"></span>
    <span class="tag" group-tag="Removal of unwanted variance"></span>
    [](./notebooks/quick_start/api_overview)
    - Go through the typical steps of an scvi-tools workflow
*   - [](./notebooks/quick_start/data_loading)
    - Load, preprocess, and register data for use with scvi-tools
*   - <span class="tag" group-tag="R"></span>
    [](./notebooks/quick_start/python_in_R)
    - Perform basic Python operations in an R environment

:::

:::

:::{dropdown} scRNA-seq

Go through the workflows of models used for analysis of scRNA-seq datasets, including scVI, scANVI, AutoZI, AmortizedLDA, LinearSCVI, and contrastiveVI

:::{list-table}
:widths: auto
:align: left

*   - <span class="tag" group-tag="Analysis"></span>
    <span class="tag" group-tag="Differential comparison"></span>
    <span class="tag" group-tag="Dimensionality reduction"></span>
    <span class="tag" group-tag="Removal of unwanted variance"></span>
    [](./notebooks/scrna/scVI_DE_worm)
    - Perform DE analysis on C. elegans data with scVI to quantify differences in gene expression between groups of cells
*   - <span class="tag" group-tag="Analysis"></span>
    <span class="tag" group-tag="Differential comparison"></span>
    <span class="tag" group-tag="Dimensionality reduction"></span>
    <span class="tag" group-tag="Removal of unwanted variance"></span>
    <span class="tag" group-tag="R"></span>
    <span class="tag" group-tag="Integration"></span>
    [](./notebooks/scrna/scvi_in_R)
    - Use basic scvi-tools functionality in R including integration of datasets
*   - <span class="tag" group-tag="Dimensionality reduction"></span>
    <span class="tag" group-tag="Removal of unwanted variance"></span>
    [](./notebooks/scrna/MrVI_tutorial)
    - <span class="tag" group-tag="Dimensionality reduction"></span>
    <span class="tag" group-tag="Removal of unwanted variance"></span>
    Analyze multi-sample scRNA-seq data with MrVI
*   - <span class="tag" group-tag="Dimensionality reduction"></span>
    <span class="tag" group-tag="Removal of unwanted variance"></span>
    [](./notebooks/scrna/MrVI_tutorial)
    - Preprocess data and train the MrVI model
*   - <span class="tag" group-tag="Dimensionality reduction"></span>
    <span class="tag" group-tag="Removal of unwanted variance"></span>
    [](./notebooks/scrna/MrVI_tutorial)
    - Visualize MrVI’s u latent space
*   - <span class="tag" group-tag="Dimensionality reduction"></span>
    <span class="tag" group-tag="Removal of unwanted variance"></span>
    <span class="tag" group-tag="Analysis"></span>
    <span class="tag" group-tag="Differential comparison"></span>
    [](./notebooks/scrna/MrVI_tutorial)
    - Perform differential expression and differential abundance analysis
*   - <span class="tag" group-tag="Analysis"></span>
    <span class="tag" group-tag="Dimensionality reduction"></span>
    <span class="tag" group-tag="Removal of unwanted variance"></span>
    <span class="tag" group-tag="Integration"></span>
    [](./notebooks/scrna/scanvi_fix)
    - Compare scANVI to other models following a bug fix in scvi-tools 1.1.0
*   - <span class="tag" group-tag="Transfer learning"></span>
    <span class="tag" group-tag="Dimensionality reduction"></span>
    <span class="tag" group-tag="Removal of unwanted variance"></span>
    [](./notebooks/scrna/seed_labeling)
    - Create seed labels and transfer cell type annotations to an entire dataset
*   - [](./notebooks/scrna/cellassign_tutorial)
    - Use CellAssign to assign cell types using only knowledge of marker genes
*   - <span class="tag" group-tag="Dimensionality reduction"></span>
    <span class="tag" group-tag="Removal of unwanted variance"></span>
    [](./notebooks/scrna/linear_decoder)
    - Fit an LDVAE model to scRNA-seq data and interpret how genes are linked to latent variables of cells
*   - [](./notebooks/scrna/AutoZI_tutorial)
    - Use the AutoZI model to enable gene-specific treatment of zero-inflation
*   - [](./notebooks/scrna/amortized_lda)
    - Run the amortized Latent Dirichlet Allocation model in scvi-tools to learn topics of an scRNA-seq dataset
*   - <span class="tag" group-tag="Dimensionality reduction"></span>
    <span class="tag" group-tag="Removal of unwanted variance"></span>
    [](./notebooks/scrna/contrastiveVI_tutorial)
    - Use contrastiveVI to isolate perturbation-induced variation in Perturb-seq data
*   - <span class="tag" group-tag="Integration"></span>
    <span class="tag" group-tag="Dimensionality reduction"></span>
    <span class="tag" group-tag="Removal of unwanted variance"></span>
    [](./notebooks/scrna/harmonization)
    - Perform integration of multiple scRNA-seq datasets both with and without cell type annotation (scVI and scANVI)
*   - <span class="tag" group-tag="Integration"></span>
    <span class="tag" group-tag="Transfer learning"></span>
    <span class="tag" group-tag="Dimensionality reduction"></span>
    <span class="tag" group-tag="Removal of unwanted variance"></span>
    [](./notebooks/scrna/tabula_muris)
    - Perform de novo integration of a labeled reference dataset with an unlabeled query dataset (label transfer)
*   - <span class="tag" group-tag="Transfer learning"></span>
    <span class="tag" group-tag="Dimensionality reduction"></span>
    <span class="tag" group-tag="Removal of unwanted variance"></span>
    [](./notebooks/scrna/scarches_scvi_tools)
    - Map cells from a query dataset to the latent space of a reference dataset with the scArches method
*   - <span class="tag" group-tag="Transfer learning"></span>
    <span class="tag" group-tag="Dimensionality reduction"></span>
    <span class="tag" group-tag="Removal of unwanted variance"></span>
    [](./notebooks/scrna/query_hlca_knn)
    - Use scANVI, scArches and scvi-hub to query the Human Lung Cell Atlas

:::
:::

:::{dropdown} scATAC-seq

Go through the workflows of models used for analysis of scATAC-seq datasets, including PeakVI, scBasset, and PoissonVI

:::{list-table}
:widths: auto
:align: left

*   - <span class="tag" group-tag="Analysis"></span>
    <span class="tag" group-tag="Differential comparison"></span>
    <span class="tag" group-tag="Dimensionality reduction"></span>
    <span class="tag" group-tag="Removal of unwanted variance"></span>
    [](./notebooks/atac/PeakVI)
    - Go through the PeakVI workflow to analyze a scATAC-seq dataset
*   - <span class="tag" group-tag="Analysis"></span>
    <span class="tag" group-tag="R"></span>
    <span class="tag" group-tag="Differential comparison"></span>
    <span class="tag" group-tag="Dimensionality reduction"></span>
    <span class="tag" group-tag="Removal of unwanted variance"></span>
    [](./notebooks/atac/peakvi_in_R)
    - Use scvi-tools functionality in R to analyze scATAC-seq data
*   - <span class="tag" group-tag="Analysis"></span>
    <span class="tag" group-tag="Dimensionality reduction"></span>
    <span class="tag" group-tag="Removal of unwanted variance"></span>
    [](./notebooks/atac/scbasset)
    - Go through the scBasset workflow to analyze a scATAC-seq dataset
*   - <span class="tag" group-tag="Integration"></span>
    <span class="tag" group-tag="Dimensionality reduction"></span>
    <span class="tag" group-tag="Removal of unwanted variance"></span>
    [](./notebooks/atac/scbasset_batch)
    - Use scBasset to integrate data across several samples
*   - <span class="tag" group-tag="Analysis"></span>
    <span class="tag" group-tag="Dimensionality reduction"></span>
    <span class="tag" group-tag="Removal of unwanted variance"></span>
    [](./notebooks/atac/PoissonVI)
    - Go through the PoissonVI workflow to analyze scATAC-seq data using quantitative fragment counts

:::
:::

:::{dropdown} scBS-seq (Single-cell bisulfite sequencing)

Use methylVI with single-cell bisulfite sequencing data

:::{list-table}
:widths: auto
:align: left

*   - <span class="tag" group-tag="Integration"></span>
    <span class="tag" group-tag="Dimensionality reduction"></span>
    <span class="tag" group-tag="Removal of unwanted variance"></span>
    [](notebooks/scbs/MethylVI_batch)
    - Correct batch effects in across different scBS-seq experiments with methylVI
*   - <span class="tag" group-tag="Integration"></span>
    <span class="tag" group-tag="Dimensionality reduction"></span>
    <span class="tag" group-tag="Removal of unwanted variance"></span>
    [](notebooks/scbs/MethylVI_batch)
    - Correct batch effects in across different scBS-seq experiments with methylVI

:::
:::

:::{dropdown} Multimodal

Use models to analyze multimodal data, including totalVI for CITE-seq analysis and MultiVI for joint analysis of scRNA-seq and scATAC-seq

:::{list-table}
:widths: auto
:align: left

*   - <span class="tag" group-tag="Analysis"></span>
    <span class="tag" group-tag="Integration"></span>
    <span class="tag" group-tag="Differential comparison"></span>
    <span class="tag" group-tag="Dimensionality reduction"></span>
    <span class="tag" group-tag="Removal of unwanted variance"></span>
    [](./notebooks/multimodal/totalVI)
    - Go through the totalVI workflow to analyze CITE-seq datasets
*   - <span class="tag" group-tag="Transfer learning"></span>
    <span class="tag" group-tag="Dimensionality reduction"></span>
    <span class="tag" group-tag="Removal of unwanted variance"></span>
    [](./notebooks/multimodal/totalVI_reference_mapping)
    - Use totalVI to train a reference model and map CITE-seq query data
*   - <span class="tag" group-tag="Integration"></span>
    <span class="tag" group-tag="Modality imputation"></span>
    <span class="tag" group-tag="Dimensionality reduction"></span>
    <span class="tag" group-tag="Removal of unwanted variance"></span>
    [](./notebooks/multimodal/cite_scrna_integration_w_totalVI)
    - Use totalVI to integrate CITE-seq and scRNA-seq datasets
*   - <span class="tag" group-tag="R"></span>
    <span class="tag" group-tag="Analysis"></span>
    <span class="tag" group-tag="Dimensionality reduction"></span>
    <span class="tag" group-tag="Removal of unwanted variance"></span>
    [](./notebooks/multimodal/totalvi_in_R)
    - Use scvi-tools functionality in R to analyze CITE-seq data
*   - <span class="tag" group-tag="Analysis"></span>
    <span class="tag" group-tag="Modality imputation"></span>
    [](./notebooks/multimodal/MultiVI_tutorial)
    - Go through the MultiVI workflow to perform joint analysis of paired and unpaired multi omic data

:::
:::

:::{dropdown} Spatial transcriptomics

Learn about models and methods to use with spatial transciptomics data, including DestVI, gimVI, Tangram, Cell2location, and Stereoscope

:::{list-table}
:widths: auto
:align: left

*   - [](./notebooks/spatial/resolVI_tutorial)
    - description
*   - <span class="tag" group-tag="Analysis"></span>
    <span class="tag" group-tag="Differential comparison"></span>
    <span class="tag" group-tag="Modality imputation"></span>
    <span class="tag" group-tag="Deconvolution"></span>
    [](./notebooks/spatial/DestVI_tutorial)
    - Perform multi-resolution analysis on spatial transcriptomics data with DestVI
*   - <span class="tag" group-tag="Analysis"></span>
    <span class="tag" group-tag="R"></span>
    <span class="tag" group-tag="Modality imputation"></span>
    <span class="tag" group-tag="Deconvolution"></span>
    [](./notebooks/spatial/DestVI_in_R)
    - Use scvi-tools functionality in R to analyze spatial transcriptomics datasets
*   - <span class="tag" group-tag="Modality imputation"></span>
    [](./notebooks/spatial/gimvi_tutorial)
    - Use gimVI to impute missing genes in spatial data
*   - <span class="tag" group-tag="Deconvolution"></span>
    [](./notebooks/spatial/tangram_scvi_tools)
    - Use Tangram to map spatial transcriptomics data
*   - <span class="tag" group-tag="Deconvolution"></span>
    [](./notebooks/spatial/stereoscope_heart_LV_tutorial)
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
*   - <span class="tag" group-tag="Analysis"></span>
    [](./notebooks/hub/scvi_hub_intro_and_download)
    - Learn how to use Hugging Face and scvi-hub to download pretrained scvi-tools models
*   - <span class="tag" group-tag="Analysis"></span>
    <span class="tag" group-tag="Differential comparison"></span>
    <span class="tag" group-tag="Dimensionality reduction"></span>
    <span class="tag" group-tag="Removal of unwanted variance"></span>
    [](./notebooks/hub/cellxgene_census_model)
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
