import json
import os

import pandas as pd
from scipy.stats import pearsonr, spearmanr
from sklearn.metrics import mean_absolute_error as mae
from sklearn.metrics import r2_score

from scvi._types import AnnOrMuData
from scvi.model._utils import REGISTRY_KEYS
from scvi.model.base import BaseModelClass

METRIC_CV_CELL = "cv_cell"
METRIC_CV_GENE = "cv_gene"
METRIC_DIFF_EXP = "diff_exp"

def _dataframe_to_markdown(df):
    # Create the header
    header = "| Index | " + " | ".join(df.columns) + " |"
    separator = "| --- | " + " | ".join("---" for _ in df.columns) + " |"

    # Format values and create rows, including the index
    rows = "\n".join(
        "| " + " | ".join(
            [str(index)] +
            [f"{value:.2f}" if not isinstance(value, int) else f"{value}" for value in row]
        ) + " |"
        for index, row in zip(df.index, df.values, strict=True)
    )

    # Combine header, separator, and rows
    return f"{header}\n{separator}\n{rows}"


def create_criticism_report(
    model: BaseModelClass,
    adata: AnnOrMuData | None = None,
    skip_metrics: list = (),
    n_samples: int = 5,
    label_key: str | None = None,
    save_folder: str | None = None,
    ) -> dict:
    """
    Helper function to compute and store criticism metrics for a model.

    Parameters
    ----------
    model
        Model to compute metrics on.
    adata
        AnnOrMuData to compute metrics on.
    skip_metrics
        List of metrics to skip. Can be one of ["cv_cell", "cv_gene", "diff_exp"].
    n_samples
        Number of samples to use for posterior predictive check.
    label_key
        Key in adata.obs to use as cell type labels. If None, will use the original label key from
        the model.
    save_folder
        Path to folder for storing the metrics. Preferred to store in save_path folder of model.
    """
    from ._ppc import PosteriorPredictiveCheck as PPC

    adata = model._validate_anndata(adata)
    models_dict = {"model": model}
    ppc = PPC(adata, models_dict, n_samples=n_samples)
    # run ppc+cv
    if METRIC_CV_CELL not in skip_metrics:
        ppc.coefficient_of_variation("features")
        md_cell_wise_cv = _cv_metrics(ppc, model, cell_wise=True)
    else:
        md_cell_wise_cv = ""
    if METRIC_CV_GENE not in skip_metrics:
        ppc.coefficient_of_variation("cells")
        md_gene_wise_cv = _cv_metrics(ppc, model, cell_wise=False)
    else:
        md_gene_wise_cv = ""
    # run diff_exp
    if METRIC_DIFF_EXP not in skip_metrics:
        labels_state_registry = model.adata_manager.get_state_registry(REGISTRY_KEYS.LABELS_KEY)
        if label_key is None and labels_state_registry.original_key!='_scvi_labels':
            label_key = labels_state_registry.original_key
        ppc.differential_expression(de_groupby=label_key, p_val_thresh=0.2)
        summary_df = ppc.metrics["diff_exp"]["summary"].set_index('group')
        summary_df = summary_df.drop(columns=['model'])
        summary_df = summary_df.sort_values(by='n_cells', ascending=False)
        md_de = _dataframe_to_markdown(summary_df)

    markdown_dict = {
        "cell_wise_cv": md_cell_wise_cv,
        "gene_wise_cv": md_gene_wise_cv,
        "diff_exp": md_de
    }

    save_path = os.path.join(save_folder, "metrics.json")

    with open(save_path, "w") as f:
        json.dump(markdown_dict, f, indent=4)


def _cv_metrics(ppc, model, cell_wise: bool = True):
    """
    Compute metrics of coefficient of variation.

    Parameters
    ----------
    model_name
        Name of the model
    cell_wise
        Whether to plot the cell-wise or gene-wise metric
    plt_type
        The type of plot to generate.
    model
        Model used for generating metric to split validation and training data
    indices
        Indices of the data to plot, if None, all data will be plotted
    """
    metric = METRIC_CV_CELL if cell_wise is True else METRIC_CV_GENE
    model_metric = ppc.metrics[metric]["model"].values
    raw_metric = ppc.metrics[metric]["Raw"].values

    # log mae, pearson corr, spearman corr, R^2
    if cell_wise:
        indices = [model.train_indices, model.validation_indices]

        mae_values = [
            mae(model_metric[indices[0]], raw_metric[indices[0]]),
            mae(model_metric[indices[1]], raw_metric[indices[1]])
        ]
        pearsonr_values = [
            pearsonr(model_metric[indices[0]], raw_metric[indices[0]])[0],
            pearsonr(model_metric[indices[1]], raw_metric[indices[1]])[0]
        ]
        spearmanr_values = [
            spearmanr(model_metric[indices[0]], raw_metric[indices[0]])[0],
            spearmanr(model_metric[indices[1]], raw_metric[indices[1]])[0]
        ]
        r2_score_values = [
            r2_score(model_metric[indices[0]], raw_metric[indices[0]]),
            r2_score(model_metric[indices[1]], raw_metric[indices[1]])
        ]
        metric_table = (
            "| Metric                  | Training Value | Validation Value |\n"
            "|-------------------------|----------------|------------------|\n"
            f"| Mean Absolute Error | {mae_values[0]:.2f}  | {mae_values[1]:.2f}           |\n"
            f"| Pearson Correlation | {pearsonr_values[0]:.2f}  | {pearsonr_values[1]:.2f}  |\n"
            f"| Spearman Correlation | {spearmanr_values[0]:.2f} | {spearmanr_values[1]:.2f}  |\n"
            f"| R² (R-Squared) | {r2_score_values[0]:.2f}  | {r2_score_values[1]:.2f}      |"
        )
    else:
        indices = [model.train_indices, model.validation_indices]

        mae_values = mae(model_metric, raw_metric)
        pearsonr_values = pearsonr(model_metric, raw_metric)[0]
        spearmanr_values = spearmanr(model_metric, raw_metric)[0]
        r2_score_values = r2_score(model_metric, raw_metric)
        metric_table = (
            "| Metric                  | Training Value |\n"
            "|-------------------------|----------------|\n"
            f"| Mean Absolute Error | {mae_values:.2f}   |\n"
            f"| Pearson Correlation | {pearsonr_values:.2f}  |\n"
            f"| Spearman Correlation | {spearmanr_values:.2f} |\n"
            f"| R² (R-Squared) | {r2_score_values:.2f}  |"
        )

    return metric_table
