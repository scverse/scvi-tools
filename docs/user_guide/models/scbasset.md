# ScBasset

**ScBasset** [^ref1] (Python class {class}`~scvi.external.SCBASSET`) posits a sequence-based method for representation learning of scATAC-seq data.

The advantages of ScBasset are:

-   Sequence representations allow for TF motif discovery and other sequence-based analyses.
-   ScBasset is fast and scalable.

The limitations of ScBasset include:

-   ScBasset cannot currently leverage unobserved data and thus cannot curently be used for tranfer learning tasks.

## Overview

TODO

[^ref1]:
    Yuan Han and David R. Kelley (2022),
    _scBasset: sequence-based modeling of single-cell ATAC-seq using convolutional neural networks_,
    [Nature Methods](https://www.nature.com/articles/s41592-022-01562-8).
