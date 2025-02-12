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
*   - <span class="tag" group-tag="Analysis"></span>
    - <span class="tag" group-tag="Integration"></span>
    - <span class="tag" group-tag="Transfer earning"></span>
    - <span class="tag" group-tag="R"></span>
:::

<!--
<span class="tag" group-tag="Analysis"></span> Analysis Tutorials
<span class="tag" group-tag="Integration"></span> Integration Tutorials
<span class="tag" group-tag="Transfer earning"></span> Transfer learning Tutorials
<span class="tag" group-tag="R"></span> R Tutorials
-->

:::{dropdown} Quick start

[](./notebooks/quick_start/api_overview.md)
[](./notebooks/quick_start/data_loading.md)
[](./notebooks/quick_start/python_in_R.md)

:::

:::{dropdown} scRNA-seq

[](./notebooks/scrna/harmonization.md)
[](./notebooks/scrna/scvi_in_R.md)
[](./notebooks/scrna/tabula_muris.md)
[](./notebooks/scrna/scarches_scvi_tools.md)
[](./notebooks/scrna/query_hlca_knn.md)
[](./notebooks/scrna/seed_labeling.md)
[](./notebooks/scrna/linear_decoder.md)
[](./notebooks/scrna/AutoZI_tutorial.md)
[](./notebooks/scrna/cellassign_tutorial.md)
[](./notebooks/scrna/amortized_lda.md)
[](./notebooks/scrna/scVI_DE_worm.md)
[](./notebooks/scrna/contrastiveVI_tutorial.md)
[](./notebooks/scrna/scanvi_fix.md)
[](./notebooks/scrna/MrVI_tutorial.md)

:::

:::{dropdown} scATAC-seq

[](./notebooks/atac/PeakVI.md)
[](./notebooks/atac/peakvi_in_R.md)
[](./notebooks/atac/scbasset.md)
[](./notebooks/atac/scbasset_batch.md)
[](./notebooks/atac/PoissonVI.md)

:::

:::{dropdown} scBS-seq (Single-cell bisulfite sequencing)

[](notebooks/scbs/MethylVI_batch.md)

:::

:::{dropdown} Multimodal

[](./notebooks/multimodal/totalVI.md)
[](./notebooks/multimodal/cite_scrna_integration_w_totalVI.md)
[](./notebooks/multimodal/totalVI_reference_mapping.md)
[](./notebooks/multimodal/totalvi_in_R.md)
[](./notebooks/multimodal/MultiVI_tutorial.md)

:::

:::{dropdown} Spatial transcriptomics

[](./notebooks/spatial/resolVI_tutorial.md)
[](./notebooks/spatial/DestVI_tutorial.md)
[](./notebooks/spatial/DestVI_in_R.md)
[](./notebooks/spatial/gimvi_tutorial.md)
[](./notebooks/spatial/tangram_scvi_tools.md)
[](./notebooks/spatial/stereoscope_heart_LV_tutorial.md)
[](./notebooks/spatial/cell2location_lymph_node_spatial_tutorial.md)

:::

:::{dropdown} Model hub

[](./notebooks/hub/cellxgene_census_model.md)
[](./notebooks/hub/scvi_hub_intro_and_download.md)
[](./notebooks/hub/scvi_hub_upload_and_large_files.md)
[](./notebooks/hub/minification.md)

:::

:::{dropdown} Hyperparameter tuning

[](./notebooks/tuning/autotune_scvi.md)
[](./notebooks/tuning/autotune_new_model.md)

:::

:::{dropdown} Development

[](./notebooks/dev/data_tutorial.md)
[](./notebooks/dev/module_user_guide.md)
[](./notebooks/dev/model_user_guide.md)

:::

[forum]: https://discourse.scverse.org/
[ecosystem]: https://scvi-tools.org/ecosystem
