# Tangram

**Tangram** {cite:p}`Biancalani21` (Python class {class}`~scvi.external.Tangram`) maps single-cell RNA-seq data to spatial data, permitting deconvolution of cell types in spatial data like Visium.

This is a reimplementation of Tangram, which can originally be found [here](https://github.com/broadinstitute/Tangram).

## Overview

Tangram learns a matrix $M$ with shape ($n_{sc} \times n_{sp}$), in which each row sums to 1. Thus this matrix can be viewed as a map from single cells to the spatial observations.

```{topic} Tutorials:

-   {doc}`/tutorials/notebooks/spatial/tangram_scvi_tools`
```
