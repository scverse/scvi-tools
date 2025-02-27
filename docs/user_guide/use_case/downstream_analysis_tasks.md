# Perform downstream analysis tasks of scvi-tools models

scvi-tools provides useful for exploring and understanding the learned latent representations, as well as for interpreting various aspects of your single-cell dataset.

1. Latent Space Exploration
Visualization: You can visualize the learned latent representation of cells to observe patterns, clusters, and structure in the data.
We usually apply Scanpy's UMAP, while this provides a visualization these methods can't fully capture the high dimensional data. See https://www.nature.com/articles/s41592-024-02301-x for a discussion of these methods.

    Clustering: After visualizing the latent space, you can perform clustering (e.g., Keiden, Kmeans, hierarchical clustering) to identify distinct groups of cells with similar gene expression profiles.

For example:

:::{figure} /\_static/img/umap.png
:align: center
:alt: General UMAP plot
:class: img-fluid
:::
2. Differential Expression (DE) Analysis
Gene Expression Comparisons: Most models allows you to perform differential expression analysis between different clusters or conditions (e.g., different cell types or experimental conditions).
You can compare the expression of genes between clusters to identify which genes are differentially expressed. To correct for covariates (e.g. sex in a comparison between healthy and diseased cells) and for a more robust approach, we recommend using pseudobulk DE methods. These models allow computation of differentially expressed genes on minified representations of the data and can leverage reference data for DE computation.
```python
differential_expression = scvi.model.SCVI().differential_expression()
```
Log-fold Change (LFC) and p-values are typically used to assess which genes have significant expression differences between groups.
3. Cell Type Identification
Mapping to Known Labels: After training a model, you can use the latent space to assign cells to known or predicted cell types. We recommend a KNN classifier a latent space. You can compare how well the learned latent space clusters cells and match them to known cell-type annotations.
If you have labeled data (e.g., cell types), you can assess how well the model’s clusters correspond to these labels.

    Marker Gene Identification: By identifying highly variable genes or differential expression, you can associate particular clusters with specific biological markers.
4. Latent Factor Analysis
Study Latent Factors: Probabilistic models learns a set of latent variables (such as cell-level and gene-level factors), which can be analyzed to understand the underlying biological drivers of gene expression. Especially helpful are linear models like LDVAE.
You can analyze latent factors to identify cell-specific effects or other latent variables that might represent biological signals such as cell cycle phases, batch effects, or experimental conditions.
5. Predictive Modeling and Imputation
Gene Expression/Protein Imputation: Models can be used to impute missing or dropout gene expression values. This is useful when dealing with sparse data, as these models can recover co-expression patterns learned from the data.
Prediction: After training, most models can also be used to predict gene expression under different conditions or experimental setups. This is called counterfactual analysis and we describe this approach in depth for MrVI.
We can get those with get_normalized_expression function and plot it, per gene over the already generated latent embedding UMAP.
6. Batch Effect Removal
Identifying and Correcting Batch Effects: Conditional models like SCVI allow you to account for and remove batch effects during training. You can use the model to test if there is any batch-related structure in the latent space or gene expression by using batch embeddings that provide low dimensional representations of batch effects.
You can evaluate how much of the variability in batch effect is due to biological factors versus technical factors like batch.
7. Trajectory Inference
Cell Fate Mapping: For datasets that span multiple conditions or time points, you can use the learned latent space to infer the trajectory of cells through developmental processes or dynamic processes like differentiation.
While our models don’t directly model dynamic trajectories, the latent representation it learns can be used as input for trajectory inference methods (like Monocle Slingshot) or optimal transport (like moscot) to infer cell differentiation paths.
8. Custom Annotations and Metadata
Incorporating Metadata: After training, models allow you to associate external metadata (such as cell type annotations or experimental conditions) with your cells and then visualize how these factors influence gene expression or the latent space.
You can annotate cells with additional biological or experimental information and analyze how they correlate with the learned latent variables.
9. Evaluation of Model Performance
Cross-validation: Perform cross-validation to check for overfitting and assess the model's generalizability on validation cells.
Comparing Models: You can compare different models trained with different configurations or hyperparameters to assess which model performs better. This is especially useful in selecting the best model for downstream tasks like differential expression or clustering.
We usually use scib-metrics to compare different models latent space and allow for posterior predictive checks within scvi-criticism to perform evaluation of the model generations.

In summary:
scvi-tools provides a broad set of downstream analysis capabilities, including differential expression analysis, cell type identification, latent factor exploration, trajectory inference, and batch effect correction, among others. By using scvi-tools’s probabilistic framework, you can explore complex patterns in single-cell RNA-seq data, visualize latent representations, impute missing data, and integrate metadata to gain deeper insights into cellular behaviors. These tools are crucial for understanding biological processes, making scvi-tools a versatile platform for single-cell genomics analysis.
