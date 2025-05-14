# scRNA-seq

```{toctree}
:maxdepth: 1

notebooks/scrna/harmonization
notebooks/scrna/MrVI_tutorial
notebooks/scrna/scVI_DE_worm
notebooks/scrna/cellassign_tutorial
notebooks/scrna/linear_decoder
notebooks/scrna/contrastiveVI_tutorial
notebooks/scrna/amortized_lda
notebooks/scrna/AutoZI_tutorial
notebooks/scrna/sysVI
notebooks/scrna/seed_labeling
notebooks/scrna/scanvi_fix
notebooks/scrna/tabula_muris
notebooks/scrna/scarches_scvi_tools
notebooks/scrna/query_hlca_knn
notebooks/scrna/decipher_tutorial
notebooks/scrna/velovi
```

```{customcard}
:path: notebooks/scrna/harmonization
:tags: Analysis, Integration, Transfer-learning, Dimensionality-reduction, Removal-of-variance

Perform integration of multiple scRNA-seq datasets both with and without cell type annotation (scVI and scANVI)
```

```{customcard}
:path: notebooks/scrna/MrVI_tutorial
:tags: Analysis, Differential-comparison, Dimensionality-reduction, Removal-of-variance

Analyze multi-sample scRNA-seq data with MrVI
```

```{customcard}
:path: notebooks/scrna/scVI_DE_worm
:tags: Analysis, Differential-comparison, Dimensionality-reduction, Removal-of-variance

Perform DE analysis on C. elegans data with scVI to quantify differences in gene expression between groups of cells
```

```{customcard}
:path: notebooks/scrna/cellassign_tutorial
:tags: Modality-imputation

Use CellAssign to assign cell types using only knowledge of marker genes
```

```{customcard}
:path: notebooks/scrna/linear_decoder
:tags: Dimensionality-reduction, Removal-of-variance, Linear-model

Fit an LDVAE model to scRNA-seq data and interpret how genes are linked to latent variables of cells
```

```{customcard}
:path: notebooks/scrna/contrastiveVI_tutorial
:tags: Analysis, Dimensionality-reduction, Removal-of-variance

Use contrastiveVI to isolate perturbation-induced variation in Perturb-seq data
```

```{customcard}
:path: notebooks/scrna/amortized_lda
:tags: Analysis, Topic-modeling

Run the amortized Latent Dirichlet Allocation model in scvi-tools to learn topics of an scRNA-seq dataset
```

```{customcard}
:path: notebooks/scrna/AutoZI_tutorial
:tags: Analysis, Dimensionality-reduction, Removal-of-variance

Use the AutoZI model to enable gene-specific treatment of zero-inflation
```

```{customcard}
:path: notebooks/scrna/sysVI
:tags: Integration, Dimensionality-reduction, Removal-of-variance

Integrate scRNA-seq datasets with substantial batch effects.
```

```{customcard}
:path: notebooks/scrna/seed_labeling
:tags: Transfer-learning, Dimensionality-reduction, Removal-of-variance

Create seed labels and transfer cell type annotations to an entire dataset
```

```{customcard}
:path: notebooks/scrna/scanvi_fix
:tags: Analysis, Integration, Transfer-learning, Dimensionality-reduction, Removal-of-variance

Compare scANVI to other models following a bug fix in scvi-tools 1.1.0
```

```{customcard}
:path: notebooks/scrna/tabula_muris
:tags: Integration, Transfer-learning, Dimensionality-reduction, Removal-of-variance

Perform de novo integration of a labeled reference dataset with an unlabeled query dataset (label transfer)
```

```{customcard}
:path: notebooks/scrna/scarches_scvi_tools
:tags: Integration, Transfer-learning, Dimensionality-reduction, Removal-of-variance

Map cells from a query dataset to the latent space of a reference dataset with the scArches method
```

```{customcard}
:path: notebooks/scrna/query_hlca_knn
:tags: Integration, Transfer-learning, Dimensionality-reduction, Removal-of-variance

Use scANVI, scArches, and scvi-hub to query the Human Lung Cell Atlas
```

```{customcard}
:path: notebooks/scrna/decipher_tutorial
:tags:

TODO: add description
```

```{customcard}
:path: notebooks/scrna/velovi
:tags:

TODO: add description
```
