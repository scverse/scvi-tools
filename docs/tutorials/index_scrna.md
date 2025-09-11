# scRNA-seq

```{toctree}
:maxdepth: 1

notebooks/scrna/harmonization
notebooks/scrna/MrVI_tutorial
notebooks/scrna/MrVI_tutorial_torch
notebooks/scrna/scanvi_fix
notebooks/scrna/seed_labeling
notebooks/scrna/tabula_muris
notebooks/scrna/scVI_DE_worm
notebooks/scrna/cellassign_tutorial
notebooks/scrna/contrastiveVI_tutorial
notebooks/scrna/linear_decoder
notebooks/scrna/amortized_lda
notebooks/scrna/AutoZI_tutorial
notebooks/scrna/sysVI
notebooks/scrna/decipher_tutorial
notebooks/scrna/velovi
notebooks/scrna/Tahoe100_mrVI
notebooks/scrna/Tahoe100_mrVI_Jax
```

```{customcard}
:path: notebooks/scrna/harmonization
:tags: Analysis, Integration

Perform integration of multiple scRNA-seq datasets both with and without cell type annotation (scVI and scANVI)
```

```{customcard}
:path: notebooks/scrna/MrVI_tutorial
:tags: Analysis, Differential-comparison, Dimensionality-reduction, Removal-of-variance

Analyze multi-sample scRNA-seq data with MrVI in Jax
```

```{customcard}
:path: notebooks/scrna/MrVI_tutorial_torch
:tags: Analysis, Differential-comparison, Dimensionality-reduction, Removal-of-variance

Analyze multi-sample scRNA-seq data with MrVI in PyTorch
```

```{customcard}
:path: notebooks/scrna/scanvi_fix
:tags: Analysis

Compare scANVI to other models following a bug fix in scvi-tools 1.1.0
```

```{customcard}
:path: notebooks/scrna/seed_labeling
:tags: Transfer-learning

Create seed labels and transfer cell type annotations to an entire dataset
```

```{customcard}
:path: notebooks/scrna/tabula_muris
:tags: Integration, Transfer-learning, Analysis

Perform de novo integration of a labeled reference dataset with an unlabeled query dataset (label transfer)
```

```{customcard}
:path: notebooks/scrna/scVI_DE_worm
:tags: Analysis, Differential-comparison, Dimensionality-reduction, Integration

Perform DE analysis on C. elegans data with scVI to quantify differences in gene expression between groups of cells
```

```{customcard}
:path: notebooks/scrna/cellassign_tutorial
:tags: Analysis

Use CellAssign to assign cell types using only knowledge of marker genes
```

```{customcard}
:path: notebooks/scrna/contrastiveVI_tutorial
:tags: Analysis, Integration, Removal-of-variance

Use contrastiveVI to isolate perturbation-induced variation in Perturb-seq data
```

```{customcard}
:path: notebooks/scrna/linear_decoder
:tags: Dimensionality-reduction, Removal-of-variance, Linear-model

Fit an LDVAE model to scRNA-seq data and interpret how genes are linked to latent variables of cells
```

```{customcard}
:path: notebooks/scrna/amortized_lda
:tags: Analysis, Topic-modeling

Run the amortized Latent Dirichlet Allocation model in scvi-tools to learn topics of an scRNA-seq dataset
```

```{customcard}
:path: notebooks/scrna/AutoZI_tutorial
:tags: Analysis

Use the AutoZI model to enable gene-specific treatment of zero-inflation
```

```{customcard}
:path: notebooks/scrna/sysVI
:tags: Integration, Analysis

Integrate scRNA-seq datasets with substantial batch effects.
```

```{customcard}
:path: notebooks/scrna/decipher_tutorial
:tags: Analysis, Integration

Use Decipher to jointly analyze samples from distinct conditions.
```

```{customcard}
:path: notebooks/scrna/velovi
:tags: Analysis

Use VeloVI to estimate RNA velocity.
```

```{customcard}
:path: notebooks/scrna/Tahoe100_mrVI
:tags: Analysis, Differential-comparison, Dimensionality-reduction, Removal-of-variance

Analyze Tahoe100M cells dataset with MrVI in PyTorch
```

```{customcard}
:path: notebooks/scrna/Tahoe100_mrVI_Jax
:tags: Analysis, Differential-comparison, Dimensionality-reduction, Removal-of-variance

Analyze Tahoe100M cells dataset with MrVI in Jax
```
