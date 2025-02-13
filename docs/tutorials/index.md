# Tutorials

The easiest way to get familiar with scvi-tools is to follow along with our tutorials.
Many are also designed to work seamlessly in Google Colab, a free cloud computing platform.
Tutorials by default work with the latest installable version of scvi-tools. To view older tutorials,
change the documentation version using the tab at the bottom of the left sidebar.

:::{note}
For questions about using scvi-tools, or broader questions about modeling data, please use our [forum]. Checkout the [ecosystem] for additional models powered by scvi-tools.
:::

```{eval-rst}
.. |integration-tag| replace:: <span class="tag" group-tag="Integration">Integration</span>
.. |analysis-tag| replace:: <span class="tag" group-tag="Analysis">Analysis</span>
.. |transfer-learning-tag| replace:: <span class="tag" group-tag="Transfer learning"></span>
.. |r-tag| replace:: <span class="tag" group-tag="R"></span>
```

:::{list-table} Tutorial Groups
:widths: auto
:align: center

*   - Analysis
    - Integration
    - Transfer learning
    - R
*   - |integration-tag|
    - |analysis-tag|
    - |transfer-learning-tag|
    - |r-tag|
:::

<!--
<span class="tag" group-tag="Analysis"></span> Analysis Tutorials
<span class="tag" group-tag="Integration"></span> Integration Tutorials
<span class="tag" group-tag="Transfer earning"></span> Transfer learning Tutorials
<span class="tag" group-tag="R"></span> R Tutorials
-->

:::{dropdown} Quick start

- [](./notebooks/quick_start/api_overview)
- [](./notebooks/quick_start/data_loading)
- |r-tag| [](./notebooks/quick_start/python_in_R)

:::

:::{dropdown} scRNA-seq

- |integration-tag| [](./notebooks/scrna/harmonization)
- |r-tag| |integration-tag| [](./notebooks/scrna/scvi_in_R)
- |integration-tag| [](./notebooks/scrna/tabula_muris)
- [](./notebooks/scrna/scarches_scvi_tools)
- [](./notebooks/scrna/query_hlca_knn)
- [](./notebooks/scrna/seed_labeling)
- [](./notebooks/scrna/linear_decoder)
- [](./notebooks/scrna/AutoZI_tutorial)
- [](./notebooks/scrna/cellassign_tutorial)
- [](./notebooks/scrna/amortized_lda)
- [](./notebooks/scrna/scVI_DE_worm)
- [](./notebooks/scrna/contrastiveVI_tutorial)
- [](./notebooks/scrna/scanvi_fix)
- [](./notebooks/scrna/MrVI_tutorial)

:::

:::{dropdown} scATAC-seq

- |analysis-tag| [](./notebooks/atac/PeakVI)
- |r-tag| {analysis-tag} [](./notebooks/atac/peakvi_in_R)
- |analysis-tag| [](./notebooks/atac/scbasset)
- [](./notebooks/atac/scbasset_batch)
- |analysis-tag| [](./notebooks/atac/PoissonVI)

:::

:::{dropdown} scBS-seq (Single-cell bisulfite sequencing)

- |integration-tag| [](notebooks/scbs/MethylVI_batch)

:::

:::{dropdown} Multimodal

- |analysis-tag| [](./notebooks/multimodal/totalVI)
- |integration-tag| [](./notebooks/multimodal/cite_scrna_integration_w_totalVI)
- [](./notebooks/multimodal/totalVI_reference_mapping)
- |analysis-tag| [](./notebooks/multimodal/totalvi_in_R)
- |analysis-tag| [](./notebooks/multimodal/MultiVI_tutorial)

:::

:::{dropdown} Spatial transcriptomics

- [](./notebooks/spatial/resolVI_tutorial)
- [](./notebooks/spatial/DestVI_tutorial)
- |r-tag| [](./notebooks/spatial/DestVI_in_R)
- [](./notebooks/spatial/gimvi_tutorial)
- [](./notebooks/spatial/tangram_scvi_tools)
- [](./notebooks/spatial/stereoscope_heart_LV_tutorial)
- [](./notebooks/spatial/cell2location_lymph_node_spatial_tutorial)

:::

:::{dropdown} Model hub

- [](./notebooks/hub/cellxgene_census_model)
- [](./notebooks/hub/scvi_hub_intro_and_download)
- [](./notebooks/hub/scvi_hub_upload_and_large_files)
- [](./notebooks/hub/minification)

:::

:::{dropdown} Hyperparameter tuning

- [](./notebooks/tuning/autotune_scvi)
- [](./notebooks/tuning/autotune_new_model)

:::

:::{dropdown} Development

- [](./notebooks/dev/data_tutorial)
- [](./notebooks/dev/module_user_guide)
- [](./notebooks/dev/model_user_guide)

:::

[forum]: https://discourse.scverse.org/
[ecosystem]: https://scvi-tools.org/ecosystem
