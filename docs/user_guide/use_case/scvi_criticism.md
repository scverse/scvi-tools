# SCVI Criticism

:::{note}
This page is under construction.
:::

SCVI-Criticism is a tool in the scvi-tools ecosystem that helps assess the performance and quality of single-cell RNA sequencing (scRNA-seq) data analysis using generative models like scVI (single-cell Variational Inference). It provides a framework for evaluating the models' predictions and helps researchers understand the strengths and weaknesses of their models. One of its main advantages is that it allows for robust evaluation using real or simulated datasets, offering insight into model robustness, overfitting, and predictive performance.

Underneath the hood it uses posterior predictive checks or PPC ({class}`scvi.criticism.PosteriorPredictiveCheck`) for comparing scRNA-seq generative models.

The method works well for any SCVI model out of the box to provide insights into the quality of predictions and model evaluation.

### There are few metrics we calculate in order to achieve that:
- **Cell Wise Coefficient of variation:** The cell-wise coefficient of variation summarizes how well variation between different cells is preserved by the generated model expression. Below a squared Pearson correlation coefficient of 0.4 , we would recommend not to use generated data for downstream analysis, while the generated latent space might still be useful for analysis.
- **Gene Wise Coefficient of variation:** The gene-wise coefficient of variation summarizes how well variation between different genes is preserved by the generated model expression. This value is usually quite high.
- **Differential expression metric:** The differential expression metric provides a summary of the differential expression analysis between cell types or input clusters. We provide the F1-score, Pearson Correlation Coefficient of Log-Foldchanges, Spearman Correlation Coefficient, and Area Under the Precision Recall Curve (AUPRC) for the differential expression analysis using Wilcoxon Rank Sum test for each cell-type.

### Example of use:
We can compute and compare 2 models PPC's by:
```python
models_dict = {"model1": model1, "model2": model2}
ppc = PPC(adata, models_dict)
```

### Creating Report:
A good practice will be to save this report together with the model, just like we do for weights and adata.
For example, for all our [SCVI-Hub](https://huggingface.co/scvi-tools) models we attached the SCVI criticism reports.

We can create the criticism report simply by:
```python
create_criticism_report(model)
```
