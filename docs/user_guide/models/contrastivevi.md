# contrastiveVI

**contrastiveVI** [^ref1] (contrastive variational inference; Python class
{class}`~scvi.external.ContrastiveVI`) is a generative model for the contrastive analysis
of scRNA-seq count data that can subsequently be used for many common downstream tasks.

Contrastive analysis requires a _target_ (e.g., treated cells) and a _background_
(e.g., control cells) dataset, and contrastiveVI is designed to isolate the variations
enriched in target cells from variations shared with background cells.

```{topic} Tutorials:

-   {doc}`/tutorials/notebooks/scrna/contrastiveVI_tutorial`
```

## Overview

:::{note}
This page is under construction.
:::

[^ref1]:
    Ethan Weinberger, Chris Lin, Su-In Lee (2023),
    _Isolating salient variations of interest in single-cell data with contrastiveVI_,
    [Nature Methods](https://www.nature.com/articles/s41592-023-01955-3).
