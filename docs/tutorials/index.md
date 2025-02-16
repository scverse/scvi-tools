# Tutorials

The easiest way to get familiar with scvi-tools is to follow along with our tutorials.
Many are also designed to work seamlessly in Google Colab, a free cloud computing platform.
Tutorials by default work with the latest installable version of scvi-tools. To view older tutorials,
change the documentation version using the tab at the bottom of the left sidebar.

:::{note}
For questions about using scvi-tools, or broader questions about modeling data, please use our [forum]. Checkout the [ecosystem] for additional models powered by scvi-tools.
:::

:::{list-table} Tutorial Groups
:widths: auto
:align: center

*   - Analysis
    - Integration
    - Transfer learning
    - R
*   - {octicon}`filled-square;1em;`
    - {octicon}`filled-square;1em;sd-color-analysis`
    - {octicon}`filled-square;1.5em;sd-color-analysis`
    -
:::

<!--
<span class="tag" group-tag="Analysis"></span> Analysis Tutorials
<span class="tag" group-tag="Integration"></span> Integration Tutorials
<span class="tag" group-tag="Transfer earning"></span> Transfer learning Tutorials
<span class="tag" group-tag="R"></span> R Tutorials
-->

:::{dropdown} Quick start

:::{list-table}
:widths: auto
:align: left

*   - [](./notebooks/quick_start/api_overview)
    - Go through the typical steps of an scvi-tools workflow
*   - [](./notebooks/quick_start/data_loading)
    - Load, preprocess, and register data for use with scvi-tools
*   - |r-tag| [](./notebooks/quick_start/python_in_R)
    - Perform basic Python operations in an R environment

:::

:::

:::{dropdown} scRNA-seq

:::{list-table}
:widths: auto
:align: left

*   - |integration-tag| [](./notebooks/scrna/harmonization)
    - Perform integration of multiple scRNA-seq datasets both with and without cell type annotation (scVI and scANVI)
*   - |r-tag| |integration-tag| [](./notebooks/scrna/scvi_in_R)
    - Use basic scvi-tools functionality in R including integration of datasets
*   - |integration-tag| [](./notebooks/scrna/tabula_muris)
    - Perform de novo integration of a labeled reference dataset with an unlabeled query dataset (label transfer)
*   - [](./notebooks/scrna/scarches_scvi_tools)
    - Map cells from a query dataset to the latent space of a reference dataset with the scArches method
*   - [](./notebooks/scrna/query_hlca_knn)
    - Use scANVI, scArches and scvi-hub to query the Human Lung Cell Atlas
*   - [](./notebooks/scrna/seed_labeling)
    - Create seed labels and transfer cell type annotations to an entire dataset
*   - [](./notebooks/scrna/linear_decoder)
    - Fit an LDVAE model to scRNA-seq data and interpret how genes are linked to latent variables of cells
*   - [](./notebooks/scrna/AutoZI_tutorial)
    - Use the AutoZI model to enable gene-specific treatment of zero-inflation
*   - [](./notebooks/scrna/cellassign_tutorial)
    - Use CellAssign to assign cell types using only knowledge of marker genes
*   - [](./notebooks/scrna/amortized_lda)
    - Run the amortized Latent Dirichlet Allocation model in scvi-tools to learn topics of an scRNA-seq dataset
*   - [](./notebooks/scrna/scVI_DE_worm)
    - Perform DE analysis on C. elegans data with scVI to quantify differences in gene expression between groups of cells
*   - [](./notebooks/scrna/contrastiveVI_tutorial)
    - Use contrastiveVI to isolate perturbation-induced variation in Perturb-seq data
*   - [](./notebooks/scrna/scanvi_fix)
    - Compare scANVI to other models following a bug fix in scvi-tools 1.1.0
*   - [](./notebooks/scrna/MrVI_tutorial)
    - Analyze multi-sample scRNA-seq data with MrVI
*   - [](./notebooks/scrna/MrVI_tutorial)
    - Preprocess data and train the MrVI model
*   - [](./notebooks/scrna/MrVI_tutorial)
    - Visualize MrVI’s u latent space
*   - [](./notebooks/scrna/MrVI_tutorial)
    - Perform differential expression and differential abundance analysis

:::

:::

:::{dropdown} scATAC-seq

:::{list-table}
:widths: auto
:align: left

*   - |analysis-tag| [](./notebooks/atac/PeakVI)
    - Go through the PeakVI workflow to analyze a scATAC-seq dataset
*   - |r-tag| {analysis-tag} [](./notebooks/atac/peakvi_in_R)
    - Use scvi-tools functionality in R to analyze scATAC-seq data
*   - |analysis-tag| [](./notebooks/atac/scbasset)
    - Go through the scBasset workflow to analyze a scATAC-seq dataset
*   - [](./notebooks/atac/scbasset_batch)
    - Use scBasset to integrate data across several samples
*   - |analysis-tag| [](./notebooks/atac/PoissonVI)
    - Go through the PoissonVI workflow to analyze scATAC-seq data using quantitative fragment counts

:::
:::

:::{dropdown} scBS-seq (Single-cell bisulfite sequencing)

:::{list-table}
:widths: auto
:align: left

*   - |integration-tag| [](notebooks/scbs/MethylVI_batch)
    - Correct batch effects in across different scBS-seq experiments with methylVI
*   - |integration-tag| [](notebooks/scbs/MethylVI_batch)
    - Correct batch effects in across different scBS-seq experiments with methylVI

:::
:::

:::{dropdown} Multimodal

:::{list-table}
:widths: auto
:align: left

*   - |analysis-tag| [](./notebooks/multimodal/totalVI)
    - Go through the totalVI workflow to analyze CITE-seq datasets
*   - |integration-tag| [](./notebooks/multimodal/cite_scrna_integration_w_totalVI)
    - Use totalVI to integrate CITE-seq and scRNA-seq datasets
*   - [](./notebooks/multimodal/totalVI_reference_mapping)
    - Use totalVI to train a reference model and map CITE-seq query data
*   - |analysis-tag| [](./notebooks/multimodal/totalvi_in_R)
    - Use scvii-tools functionality in R to analyze CITE-seq data
*   - |analysis-tag| [](./notebooks/multimodal/MultiVI_tutorial)
    - Go through the MultiVI workflow to perform joint analysis of paired and unpaired multi omic data

:::
:::

:::{dropdown} Spatial transcriptomics

:::{list-table}
:widths: auto
:align: left

*   - [](./notebooks/spatial/resolVI_tutorial)
    - description
*   - [](./notebooks/spatial/DestVI_tutorial)
    - Perform multi-resolution analysis on spatial transcriptomics data with DestVI
*   - |r-tag| [](./notebooks/spatial/DestVI_in_R)
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

:::{list-table}
:widths: auto
:align: left

*   - [](./notebooks/hub/cellxgene_census_model)
    - Perform analysis of a CELLxGENE dataset using a pretrained model from scVI-hub
*   - [](./notebooks/hub/scvi_hub_intro_and_download)
    - Learn how to use Hugging Face and scvi-hub to download pretrained scvi-tools models
*   - [](./notebooks/hub/scvi_hub_upload_and_large_files)
    - Learn how to upload pretrained scvi-tools models to Hugging Face
*   - [](./notebooks/hub/minification)
    - Minify a dataset by replacing count data with the model’s estimated parameters of the latent posterior distribution


:::
:::

:::{dropdown} Hyperparameter tuning

:::{list-table}
:widths: auto
:align: left

*   - [](./notebooks/tuning/autotune_scvi)
    - Use scvi’s autotune module to automatically find a good set of model hyperparameters
*   - [](./notebooks/tuning/autotune_new_model)
    - description

:::
:::

:::{dropdown} Development

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
