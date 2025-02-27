# scBasset

**scBasset** [^ref1] (Python class {class}`~scvi.external.SCBASSET`) posits a sequence-based method for representation learning of scATAC-seq data.

The advantages of ScBasset are:

-   Sequence representations allow for TF motif discovery and other sequence-based analyses.
-   scBasset is fast and scalable.

The limitations of scBasset include:

-   scBasset cannot currently leverage unobserved data and thus cannot curently be used for tranfer learning tasks.

## Overview

:::{note}
This page is under construction.
:::

```{topic} Tutorials:

-   {doc}`/tutorials/notebooks/atac/scbasset`
-   {doc}`/tutorials/notebooks/atac/scbasset_batch`
```

[^ref1]:
    Yuan Han and David R. Kelley (2022),
    _scBasset: sequence-based modeling of single-cell ATAC-seq using convolutional neural networks_,
    [Nature Methods](https://www.nature.com/articles/s41592-022-01562-8).
