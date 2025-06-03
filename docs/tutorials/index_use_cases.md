# Common Modelling Use Cases

```{toctree}
:maxdepth: 1

notebooks/use_cases/preprocessing
notebooks/use_cases/autotune_scvi
notebooks/use_cases/minification
notebooks/use_cases/interpretability
notebooks/use_cases/custom_dl/tiledb
notebooks/use_cases/custom_dl/lamin
notebooks/use_cases/custom_dl/ann_collection
notebooks/use_cases/multiGPU
```

```{customcard}
:path: notebooks/use_cases/preprocessing
:tags:

Learn how to preprocess various types of data for use with scvi-tools models.
```

```{customcard}
:path: notebooks/use_cases/autotune_scvi
:tags: Analysis, Integration, Dev

Automatically find optimal set of hyperparameters using autotune.
```

```{customcard}
:path: notebooks/use_cases/minification
:tags: Minification, Dev

Minify a dataset by replacing count data with the model’s estimated parameters of the latent posterior distribution
```

```{customcard}
:path: notebooks/use_cases/interpretability
:tags: Analysis

Use integrated gradient or SHAP values for model explainability
```

```{customcard}
:path: notebooks/use_cases/custom_dl/tiledb
:tags: Analysis, Custom-Dataloaders, Integration, Dev

Learn a scalable approach using TileDBDataModule dataloader to training an scVI model on Census data.
```

```{customcard}
:path: notebooks/use_cases/custom_dl/lamin
:tags: Analysis, Custom-Dataloaders, Integration, Dev

Use the Lamin MappedCollectionDataModule for a scalable approach to training an scVI model on multiple adata's.
```

```{customcard}
:path: notebooks/use_cases/custom_dl/ann_collection
:tags: Analysis, Custom-Dataloaders, Integration, Dev

Use the AnnCollection dataloader for a scalable approach to training an scVI model on multiple adata's.
```

```{customcard}
:path: notebooks/use_cases/multiGPU
:tags: Analysis

Example of how to train an SCVI model using multi GPU settings
```
