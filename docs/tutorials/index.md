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
*   - {octicon}`filled-square;1em;custom-analysis`
    - 
    - 
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
    - description
*   - [](./notebooks/quick_start/data_loading)
    - description
*   - |r-tag| [](./notebooks/quick_start/python_in_R)
    - description

:::

:::

:::{dropdown} scRNA-seq

:::{list-table}
:widths: auto
:align: left

*   - |integration-tag| [](./notebooks/scrna/harmonization)
    - description
*   - |r-tag| |integration-tag| [](./notebooks/scrna/scvi_in_R)
    - description
*   - |integration-tag| [](./notebooks/scrna/tabula_muris)
    - description
*   - [](./notebooks/scrna/scarches_scvi_tools)
    - description
*   - [](./notebooks/scrna/query_hlca_knn)
    - description
*   - [](./notebooks/scrna/seed_labeling)
    - description
*   - [](./notebooks/scrna/linear_decoder)
    - description
*   - [](./notebooks/scrna/AutoZI_tutorial)
    - description
*   - [](./notebooks/scrna/cellassign_tutorial)
    - description
*   - [](./notebooks/scrna/amortized_lda)
    - description
*   - [](./notebooks/scrna/scVI_DE_worm)
    - description
*   - [](./notebooks/scrna/contrastiveVI_tutorial)
    - description
*   - [](./notebooks/scrna/scanvi_fix)
    - description
*   - [](./notebooks/scrna/MrVI_tutorial)
    - description

:::

:::

:::{dropdown} scATAC-seq

:::{list-table}
:widths: auto
:align: left

*   - |analysis-tag| [](./notebooks/atac/PeakVI)
    - description
*   - |r-tag| {analysis-tag} [](./notebooks/atac/peakvi_in_R)
    - description
*   - |analysis-tag| [](./notebooks/atac/scbasset)
    - description
*   - [](./notebooks/atac/scbasset_batch)
    - description
*   - |analysis-tag| [](./notebooks/atac/PoissonVI)
    - description

:::
:::

:::{dropdown} scBS-seq (Single-cell bisulfite sequencing)

:::{list-table}
:widths: auto
:align: left

*   - |integration-tag| [](notebooks/scbs/MethylVI_batch)
    - description

:::
:::

:::{dropdown} Multimodal

:::{list-table}
:widths: auto
:align: left

*   - |analysis-tag| [](./notebooks/multimodal/totalVI)
    - description
*   - |integration-tag| [](./notebooks/multimodal/cite_scrna_integration_w_totalVI)
    - description
*   - [](./notebooks/multimodal/totalVI_reference_mapping)
    - description
*   - |analysis-tag| [](./notebooks/multimodal/totalvi_in_R)
    - description
*   - |analysis-tag| [](./notebooks/multimodal/MultiVI_tutorial)
    - description

:::
:::

:::{dropdown} Spatial transcriptomics

:::{list-table}
:widths: auto
:align: left

*   - [](./notebooks/spatial/resolVI_tutorial)
    - description
*   - [](./notebooks/spatial/DestVI_tutorial)
    - description
*   - |r-tag| [](./notebooks/spatial/DestVI_in_R)
    - description
*   - [](./notebooks/spatial/gimvi_tutorial)
    - description
*   - [](./notebooks/spatial/tangram_scvi_tools)
    - description
*   - [](./notebooks/spatial/stereoscope_heart_LV_tutorial)
    - description
*   - [](./notebooks/spatial/cell2location_lymph_node_spatial_tutorial)
    - description

:::
:::

:::{dropdown} Model hub

:::{list-table}
:widths: auto
:align: left

*   - [](./notebooks/hub/cellxgene_census_model)
    - description
*   - [](./notebooks/hub/scvi_hub_intro_and_download)
    - description
*   - [](./notebooks/hub/scvi_hub_upload_and_large_files)
    - description
*   - [](./notebooks/hub/minification)
    - description

:::
:::

:::{dropdown} Hyperparameter tuning

:::{list-table}
:widths: auto
:align: left

*   - [](./notebooks/tuning/autotune_scvi)
    - description
*   - [](./notebooks/tuning/autotune_new_model)
    - description

:::
:::

:::{dropdown} Development

:::{list-table}
:widths: auto
:align: left

*   - [](./notebooks/dev/data_tutorial)
    - description
*   - [](./notebooks/dev/module_user_guide)
    - description
*   - [](./notebooks/dev/model_user_guide)
    - description

:::
:::

[forum]: https://discourse.scverse.org/
[ecosystem]: https://scvi-tools.org/ecosystem
