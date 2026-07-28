# scBasset

**scBasset** {cite:p}`Yuan2022` (Python class {class}`~scvi.external.SCBASSET`) is a
sequence-based model for representation learning of scATAC-seq data. It uses the DNA
sequence of each accessible region to learn region embeddings and jointly learns cell
embeddings that reconstruct a binary accessibility matrix.

:::{warning}
SCBASSET's development is still in progress. The current scvi-tools implementation may
not fully reproduce the original implementation's results.
:::

The advantages of scBasset are:

-   Sequence representations allow for TF motif discovery and other sequence-based analyses.
-   The learned cell embeddings can be used for visualization, clustering, and batch
    integration of scATAC-seq data.
-   The model can score transcription factor activity with a motif injection procedure.

The limitations of scBasset include:

-   It expects binary accessibility data and DNA sequence encodings for the genomic
    regions.
-   The current implementation assumes fixed-length sequence inputs, following the
    original 1344 bp scBasset setting.
-   scBasset cannot currently leverage unobserved data and thus cannot currently be used
    for transfer learning tasks.
-   The built-in motif library download currently supports the human motif library used
    by the scBasset paper.

```{topic} Tutorials:

-   {doc}`/tutorials/notebooks/atac/scbasset`
-   {doc}`/tutorials/notebooks/atac/scbasset_batch`
```

## Preliminaries

scBasset uses a region-by-cell AnnData object. In a standard scATAC-seq AnnData object,
cells are observations and regions are variables, so the data are typically transposed
before setup:

```
>>> bdata = adata.transpose()
>>> SCBASSET.setup_anndata(bdata, layer="binary", dna_code_key="dna_code")
```

The registered matrix should contain binary accessibility values. The `dna_code_key`
argument points to integer-encoded DNA sequences for each region. In the transposed
object, these encodings are stored in `bdata.obsm`, one row per region. If a `batch_key`
is supplied, it is read from `bdata.var` because cells are variables in this layout.

The tutorial demonstrates creating the required sequence fields with
{func}`~scvi.data.add_dna_sequence`, which stores both raw sequence strings and integer
codes.

## Overview

scBasset is not a variational autoencoder. It is a neural network that predicts
cell-by-region accessibility from genomic sequence.

The model first converts each DNA sequence into a one-hot representation. A convolutional
neural network processes the sequence with stochastic reverse-complement augmentation,
stochastic shifts, a stem convolution, a convolutional tower, and a bottleneck dense
layer. The output is a low-dimensional embedding for each genomic region.

The model also learns:

-   a cell embedding matrix,
-   a cell-specific bias term, and
-   when a batch key is registered, a batch embedding for each cell's batch.

The accessibility logits are computed as the matrix product between region embeddings
and cell embeddings, plus the cell bias. When batches are registered, the batch embedding
is added to the cell embedding before this product.

## Inference

scBasset is trained by minimizing binary cross-entropy between predicted accessibility
logits and the observed binary accessibility matrix. The implementation also reports
AUROC during training and can add L2 regularization to the cell embedding matrix with
`l2_reg_cell_embedding`, which is useful in the batch-integration tutorial.

Training mini-batches are over regions, not cells. This follows from the region-by-cell
input layout and the sequence encoder, which processes a batch of region sequences at a
time.

## Tasks

Here we provide an overview of common tasks. Please see {class}`~scvi.external.SCBASSET`
for the full API reference.

### Cell Representation

The learned cell embedding is returned by
{meth}`~scvi.external.SCBASSET.get_latent_representation`:

```
>>> adata.obsm["X_scbasset"] = model.get_latent_representation()
```

This representation can be used for nearest-neighbor graph construction, visualization,
clustering, or integration diagnostics.

### Cell Bias

{meth}`~scvi.external.SCBASSET.get_cell_bias` returns the learned cell-specific bias
term, which reflects cell-level accessibility propensity in the reconstruction model.

### Transcription Factor Activity

{meth}`~scvi.external.SCBASSET.get_tf_activity` estimates transcription factor activity
with motif injection. The method compares model-predicted accessibility for sequences
with a known motif inserted against dinucleotide-shuffled background sequences, then
returns a cell-level activity score for the requested transcription factor.
