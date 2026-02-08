scvi_pretext = """
ScVI is a variational inference model for single-cell RNA-seq data that can learn an underlying
latent space, integrate technical batches and impute dropouts.
The learned low-dimensional latent representation of the data can be used for visualization and
clustering.

scVI takes as input a scRNA-seq gene expression matrix with cells and genes.
We provide an extensive [user guide](https://docs.scvi-tools.org/en/stable/user_guide/models/scvi.html).

- See our original manuscript for further details of the model:
[scVI manuscript](https://www.nature.com/articles/s41592-018-0229-2).
- See our manuscript on [scvi-hub](https://www.biorxiv.org/content/10.1101/2024.03.01.582887v2) how
to leverage pre-trained models.

This model can be used for fine tuning on new data using our Arches framework:
[Arches tutorial](https://docs.scvi-tools.org/en/stable/tutorials/notebooks/scrna/scarches_scvi_tools.html).
"""

scanvi_pretext = """
ScANVI is a variational inference model for single-cell RNA-seq data that can learn an underlying
latent space, integrate technical batches and impute dropouts.
In addition, to scVI, ScANVI is a semi-supervised model that can leverage labeled data to learn a
cell-type classifier in the latent space and afterward predict cell types of new data.
The learned low-dimensional latent representation of the data can be used for visualization and
clustering.

scANVI takes as input a scRNA-seq gene expression matrix with cells and genes as well as a
cell-type annotation for a subset of cells.
We provide an extensive [user guide](https://docs.scvi-tools.org/en/stable/user_guide/models/scanvi.html).

- See our original manuscript for further details of the model:
[scANVI manuscript](https://www.embopress.org/doi/full/10.15252/msb.20209620).
- See our manuscript on [scvi-hub](https://www.biorxiv.org/content/10.1101/2024.03.01.582887v2)
how to leverage pre-trained models.

This model can be used for fine tuning on new data using our Arches framework:
[Arches tutorial](https://docs.scvi-tools.org/en/stable/tutorials/notebooks/scrna/scarches_scvi_tools.html).
"""

condscvi_pretext = """
CondSCVI is a variational inference model for single-cell RNA-seq data that can learn an underlying
latent space. The predictions of the model are meant to be afterward
used for deconvolution of a second spatial transcriptomics dataset in DestVI. DestVI predicts the
cell-type proportions as well as cell type-specific activation state
in the spatial data.

CondSCVI takes as input a scRNA-seq gene expression matrix with cells and genes as well as a
cell-type annotation for all cells.
We provide an extensive [user guide](https://docs.scvi-tools.org/en/stable/user_guide/models/destvi.html)
for DestVI including a description of CondSCVI.

- See our original manuscript for further details of the model:
[DestVI manuscript](https://www.nature.com/articles/s41587-022-01272-8).
- See our manuscript on [scvi-hub](https://www.biorxiv.org/content/10.1101/2024.03.01.582887v2)
how to leverage pre-trained models.
"""

stereoscope_pretext = """
Stereoscope is a variational inference model for single-cell RNA-seq data that can learn a
cell-type specific rate of gene expression. The predictions of the model are meant to be afterward
used for deconvolution of a second spatial transcriptomics dataset in Stereoscope. Stereoscope
predicts the cell-type proportions in the spatial data.

Stereoscope takes as input a scRNA-seq gene expression matrix with cells and genes as well as a
cell-type annotation for all cells.
We provide an extensive for DestVI including a description of CondSCVI
[user guide](https://docs.scvi-tools.org/en/stable/user_guide/models/destvi.html).

- See our original manuscript for further details of the model:
[Stereoscope manuscript](https://www.nature.com/articles/s42003-020-01247-y) as well as the
[scvi-tools manuscript](https://www.nature.com/articles/s41587-021-01206-w) about implementation
details.
- See our manuscript on [scvi-hub](https://www.biorxiv.org/content/10.1101/2024.03.01.582887v2)
how to leverage pre-trained models.
"""

totalvi_pretext = """
TotalVI is a variational inference model for single-cell RNA-seq as well as protein data that can
learn an underlying latent space, integrate technical batches, impute dropouts,
and predict protein expression given gene expression or missing protein data given gene expression
and protein data for a subset of proteins.
The learned low-dimensional latent representation of the data can be used for visualization and
clustering.

TotalVI takes as input a scRNA-seq gene expression and protein expression matrix with cells and
genes.
We provide an extensive [user guide](https://docs.scvi-tools.org/en/stable/user_guide/models/totalvi.html).

- See our original manuscript for further details of the model:
[TotalVI manuscript](https://www.nature.com/articles/s41592-020-01050-x).
- See our manuscript on [scvi-hub](https://www.biorxiv.org/content/10.1101/2024.03.01.582887v2)
how to leverage pre-trained models.

This model can be used for fine tuning on new data using our Arches framework:
[Arches tutorial](https://docs.scvi-tools.org/en/stable/tutorials/notebooks/scrna/scarches_scvi_tools.html).
"""


template = """\
---
{card_data}
---

{pretext}

# Model Description

{description}

# Metrics

We provide here key performance metrics for the uploaded model, if provided by the data uploader.

<details>
<summary><strong>Coefficient of variation</strong></summary>

The cell-wise coefficient of variation summarizes how well variation between different cells is
preserved by the generated model expression. Below a squared Pearson correlation coefficient of 0.4
, we would recommend not to use generated data for downstream analysis, while the generated latent
space might still be useful for analysis.

**Cell-wise Coefficient of Variation**:

{cell_wise_cv}

The gene-wise coefficient of variation summarizes how well variation between different genes is
preserved by the generated model expression. This value is usually quite high.

**Gene-wise Coefficient of Variation**:

{gene_wise_cv}

</details>

<details>
<summary><strong>Differential expression metric</strong></summary>

The differential expression metric provides a summary of the differential expression analysis
between cell types or input clusters. We provide here the F1-score, Pearson Correlation
Coefficient of Log-Foldchanges, Spearman Correlation Coefficient, and Area Under the Precision
Recall Curve (AUPRC) for the differential expression analysis using Wilcoxon Rank Sum test for each
cell-type.

**Differential expression**:

{de_metrics}

</details>

# Model Properties

We provide here key parameters used to setup and train the model.

<details>
<summary><strong>Model Parameters</strong></summary>

These provide the settings to setup the original model:
```json
{model_init_params}
```

</details>

<details>
<summary><strong>Setup Data Arguments</strong></summary>

Arguments passed to setup_anndata of the original model:
```json
{model_setup_anndata_args}
```

</details>

<details>
<summary><strong>Data Registry</strong></summary>

Registry elements for AnnData manager:
{model_data_registry}

- **Data is Minified**: {data_is_minified}

</details>

<details>
<summary><strong>Summary Statistics</strong></summary>

{model_summary_stats}

</details>


<details>
<summary><strong>Training</strong></summary>

<!-- If your model is not uploaded with any data (e.g., minified data) on the Model Hub, then make
sure to provide this field if you want users to be able to access your training data. See the
scvi-tools documentation for details. -->
**Training data url**: {training_data_url}

If provided by the original uploader, for those interested in understanding or replicating the
training process, the code is available at the link below.

**Training Code URL**: {training_code_url}

</details>
\n
# References

{references}
"""
