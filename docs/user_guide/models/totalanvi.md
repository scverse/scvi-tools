# TotalANVI

:::{note}
This page is under construction.
:::

**TotalANVI** [^ref1] (Python class {class}`~scvi.external.TOTALANVI`) is a semi-supervised generative model of CITE-seq RNA and protein data.
Similar to how scANVI extends scVI, TotalANVI can be treated as an extension of TotalVI that can leverage cell type annotations
for a subset of the cells present in the data sets to infer the states of the rest of the cells as well as impute missing proteins expression

The advantages of TotalANVI are:

-   Comprehensive in capabilities.
-   Scalable to very large datasets (>1 million cells).

The limitations of TotalANVI include:

-   Effectively requires a GPU for fast inference.
-   May not scale to very large number of cell types.

```{topic} Tutorials:

-   Work in progress.
```

## Preliminaries

## Generative process

## Inference

## Training details

## Tasks

### Cell type label prediction


[^ref1]:
