# Perform downstream analysis tasks of SCVI models

:::{note}
In order to run scvi-tools with scanpy support, use: pip install scvi-tools[scanpy]
:::

SCVI provides useful tools for exploring and understanding the learned latent representations, as well as for interpreting various aspects of your single-cell dataset.

1. Latent Space Exploration
Visualization: You can visualize the learned latent representation of cells to observe patterns, clusters, and structure in the data.
We usually apply [scanpy's UMAP](https://scanpy.readthedocs.io/en/1.10.x/tutorials/plotting/core.html); These dimensionality reduction techniques can be applied to the model’s latent space to visualize how cells are grouped in 2D or 3D.

    Clustering: After visualizing the latent space, you can perform clustering (e.g., leiden, k-means, hierarchical clustering) to identify distinct groups of cells with similar gene expression profiles.

For example:

:::{figure} /\_static/img/umap.png
:align: center
:alt: General UMAP plot
:class: img-fluid
:::
2. Differential Expression (DE) Analysis
Gene Expression Comparisons: SCVI allows you to perform differential expression analysis between different clusters or conditions (e.g., different cell types or experimental conditions).
You can compare the expression of genes between clusters to identify which genes are differentially expressed. See more information [here](https://decoupler-py.readthedocs.io/en/latest/notebooks/bulk.html#Differential-expression-analysis)
```python
differential_expression = scvi.model.SCVI().differential_expression()
```
Log-fold Change (LFC) and p-values are typically used to assess which genes have significant expression differences between groups.
Refer to [SCVI-Hub](https://huggingface.co/scvi-tools) for use cases of DE.
3. Cell Type Identification
Mapping to Known Labels: After training a model with SCVI, you can use the latent space to assign cells to known or predicted cell types. You can compare how well SCVI clusters cells by their latent representations and match them to known biological annotations.
If you have labeled data (e.g., cell types), you can assess how well the model’s clusters correspond to these labels.

    Marker Gene Identification: By identifying highly variable genes or differential expression, you can associate particular clusters with specific biological markers.
4. Latent Factor Analysis
Study Latent Factors: SCVI learns a set of latent variables (such as cell-level and gene-level factors), which can be analyzed to understand the underlying biological drivers of gene expression.
You can analyze latent factors to identify cell-specific effects or other latent variables that might represent biological signals such as cell cycle phases, batch effects, or experimental conditions.
See {class}`scvi.module.LDVAE`
5. Predictive Modeling and Imputation
Gene Expression/Protein Imputation: SCVI can be used to impute missing or dropout gene expression values. This is useful when dealing with sparse data, as SCVI can recover hidden information based on the patterns it learns in the data.
Prediction: After training, SCVI can also be used to predict gene expression under different conditions or experimental setups.
We can get those with get_normalized_expression function, exists for most models, and plot it, per gene over the already generated latent embedding umap.
6. Batch Effect Removal
Identifying and Correcting Batch Effects: SCVI allows you to account for and remove batch effects during training. You can use the model to test if there is any batch-related structure in the latent space or gene expression.
You can evaluate how much of the variability in gene expression is due to biological factors versus technical factors like batch.
7. Trajectory Inference
Cell Fate Mapping: For datasets that span multiple conditions or time points, you can use the learned latent space to infer the trajectory of cells through developmental processes or dynamic processes like differentiation.
While SCVI doesn’t directly model dynamic trajectories, the latent representation it learns can be used as input for trajectory inference methods (like Monocle or Slingshot) to infer cell differentiation paths.
8. Custom Annotations and Metadata
Incorporating Metadata: After training, SCVI allows you to associate external metadata (such as cell type annotations or experimental conditions) with your cells and then visualize how these factors influence gene expression or the latent space.
You can annotate cells with additional biological or experimental information and analyze how they correlate with the learned latent variables.
9. Evaluation of Model Performance
Log-Likelihood: SCVI models the data through a probabilistic framework, so you can evaluate the model by looking at its log-likelihood to assess how well it fits the data.
Cross-validation: Perform cross-validation to check for overfitting and assess the model's generalizability.
10. Model Comparison
Comparing SCVI Models: You can compare different SCVI models trained with different configurations or hyperparameters to assess which model performs better. This is especially useful in selecting the best model for downstream tasks like differential expression or clustering.
We usually use to run [scib-metrics](https://github.com/YosefLab/scib-metrics) to compare different models bio conservation and batch correction metrics.

In summary:
SCVI provides a broad set of downstream analysis capabilities, including differential expression analysis, cell type identification, latent factor exploration, trajectory inference, and batch effect correction, among others. By using SCVI’s probabilistic framework, you can explore complex patterns in single-cell RNA-seq data, visualize latent representations, impute missing data, and integrate metadata to gain deeper insights into cellular behaviors. These tools are crucial for understanding biological processes, making SCVI a versatile tool for single-cell genomics analysis.
