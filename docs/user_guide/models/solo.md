# Solo

**Solo** [^ref1] (Python class {class}`~scvi.external.SOLO`) posits a flexible generative model of scRNA-seq count data that can subsequently
be used for many common downstream tasks.

The advantages of Solo are:

-   Can perform doublet detection on pre-trained {class}`~scvi.model.SCVI` models
-   Scalable to very large datasets (>1 million cells).

The limitations of Solo include:

-   For an analysis seeking to only do doublet detection, Solo will be slower than other methods.

## Overview

Solo starts with a trained {class}`~scvi.model.SCVI` instance. First Solo, simulates doublets using
the original data and second Solo trains a classifer on the model latent space.

### Doublet simulation

A simulated doublet $d_n$ is generated via the following process:

```{math}
:nowrap: true

\begin{align}
 d_n  = x_{1} + x_{2},
 \end{align}
```

where $x_{1}$ and $x_{2}$ are drawn i.i.d from the
empirical data distribution $p_{\textrm{data}}(x)$ over single-cell
transcriptomes (count data).

The number of doublets to generate is controlled by the `doublet_ratio` parameter of
{func}`~scvi.external.SOLO.from_scvi_model`.

### Classifier training

After doublet simulation, the doublets are encoded through the scVI encoder, which outputs latent
representations $z'_{1:D}$ if there are $D$ doublets.

These vectors are assigned a label of 1, while the latent representations of the original data $z_{1:N}$ are
assigned a label of 0. A simple multilayer perceptron classifier ({class}`scvi.module.Classifier`) is trained
and the doublet score for each originally observed cell is the doublet probability according to this classifier.

[^ref1]:
    Nicholas J. Bernstein, , Nicole L. Fong, Irene Lam, Margaret A. Roy, David G. Hendrickson, and David R. Kelley (2020),
    _Solo: doublet identification in single-cell RNA-Seq via semi-supervised deep learning_,
    [Cell Systems](https://www.sciencedirect.com/science/article/pii/S2405471220301952).
