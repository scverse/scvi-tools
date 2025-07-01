# CytoVI

**CytoVI** (Python class {class}`~scvi.external.CYTOVI`) is a generative model for cytometry data that leverages deep probabilistic latent variable modeling to enable denoising, imputation, integration, and differential analysis across technologies and batches.

The advantages of CytoVI are:

- Provides a denoised and batch-corrected low-dimensional cell state representations across antibody-based single cell technologies.
- Facilitates the integration of data from different antibody panels and the imputation of non-overlapping markers.
- Enables cross-technology integration (e.g. flow cytometry, mass cytometry and CITE-seq) and modality imputation, facilitating unified analysis across cytometry and transcriptomics.
- Allows for the automated identification of disease-associated cell states (label-free differential abundance analysis) and differential protein expression between groups of cells.
- Scalable to very large datasets (>20 million cells).

The limitations of CytoVI include:

- Requires at least partial feature overlap across datasets for effective integration.
- Effectively requires a GPU for fast inference.
- Assumes measurements have been corrected for fluorescent spillover and preprocessed according to standard practice.

```{topic} Tutorials:

- {doc}`/tutorials/notebooks/cytometry/CytoVI_tutorial` # note: to populate

## Preliminaries


## Inference


## Tasks
