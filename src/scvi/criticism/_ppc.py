from __future__ import annotations

import json
import warnings
from typing import TYPE_CHECKING, Literal

import numpy as np
import pandas as pd
from anndata import AnnData
from mudata import MuData
from scipy.sparse import issparse
from scipy.stats import pearsonr, spearmanr
from sklearn.metrics import (
    average_precision_score,
    precision_recall_fscore_support,
    roc_auc_score,
)
from sparse import GCXS, SparseArray
from xarray import DataArray, Dataset

from scvi.utils import dependencies

from ._constants import (
    DATA_VAR_RAW,
    DEFAULT_DE_N_TOP_GENES_OVERLAP,
    DEFAULT_DE_P_VAL_THRESHOLD,
    METRIC_CALIBRATION,
    METRIC_CV_CELL,
    METRIC_CV_GENE,
    METRIC_DIFF_EXP,
    METRIC_ZERO_FRACTION,
    UNS_NAME_RGG_PPC,
    UNS_NAME_RGG_RAW,
)

if TYPE_CHECKING:
    from scvi._types import AnnOrMuData
    from scvi.model.base import BaseModelClass

Dims = Literal["cells", "features"]


def _make_dataset_dense(dataset: Dataset) -> Dataset:
    """Make a dataset dense, converting sparse arrays to dense arrays."""
    dataset = dataset.map(lambda x: x.data.todense() if isinstance(x.data, SparseArray) else x)
    return dataset


class PosteriorPredictiveCheck:
    """
    Posterior predictive checks for comparing scRNA-seq generative models.

    Parameters
    ----------
    adata
        :class:`~AnnOrMudata` object with raw counts in either ``adata.X`` or ``adata.layers``.
    models_dict
        Dictionary of models to compare.
    count_layer_key
        Key in ``adata.layers`` to use as raw counts. If ``None``, defaults to ``adata.X``.
    n_samples
        Number of posterior predictive samples to generate.
    indices
        Indices of observations in ``adata`` to subset to before generating posterior predictive
        samples and computing metrics. If ``None``, defaults to all observations in ``adata``.
    modality
        Modality to use for posterior predictive samples. Needs to be defined if using MuData
    """

    def __init__(
        self,
        adata: AnnOrMuData,
        models_dict: dict[str, BaseModelClass],
        count_layer_key: str | None = None,
        n_samples: int = 10,
        indices: list | None = None,
        modality: str | None = None,
    ):
        if indices is not None:
            adata = adata[indices]
        self.count_layer_key = count_layer_key
        self.modality = modality
        if isinstance(adata, MuData):
            assert modality is not None, "Modality must be defined for MuData."
            self.adata = adata[modality]
            raw_counts = (
                self.adata.layers[count_layer_key] if count_layer_key is not None else self.adata.X
            )
        else:
            self.adata = adata
            raw_counts = adata.layers[count_layer_key] if count_layer_key is not None else adata.X
        # Compressed axis is rows, like csr
        if isinstance(raw_counts, np.ndarray):
            self.raw_counts = GCXS.from_numpy(raw_counts, compressed_axes=(0,))
        elif issparse(raw_counts):
            self.raw_counts = GCXS.from_scipy_sparse(raw_counts).change_compressed_axes((0,))
        else:
            raise ValueError("raw_counts must be a numpy array or scipy sparse matrix")
        self.samples_dataset = None
        self.n_samples = n_samples
        self.models = models_dict
        self.metrics = {}

        self._store_posterior_predictive_samples(indices=indices)

    def __repr__(self) -> str:
        return (
            f"--- Posterior Predictive Checks ---\n"
            f"n_samples = {self.n_samples}\n"
            f"raw_counts shape = {self.raw_counts.shape}\n"
            f"models: {list(self.models.keys())}\n"
            f"metrics: \n{self._metrics_repr()}"
        )

    def _metrics_repr(self) -> str:
        def custom_handle_unserializable(o):
            if isinstance(o, AnnData):
                return f"AnnData object with n_obs={o.n_obs}, n_vars={o.n_vars}"
            elif isinstance(o, pd.DataFrame):
                s = f"Pandas DataFrame with shape={o.shape}, "
                n_cols = 5
                if len(o.columns) > n_cols:
                    return s + f"first {n_cols} columns={o.columns[:n_cols].to_list()}"
                return s + f"columns={o.columns.to_list()}"
            elif isinstance(o, pd.Series):
                return f"Pandas Series with n_rows={len(o)}"
            return f"ERROR unserializable type: {type(o)}"

        return json.dumps(self.metrics, indent=4, default=custom_handle_unserializable)

    def _store_posterior_predictive_samples(
        self,
        batch_size: int = 32,
        indices: list[int] | None = None,
    ):
        """
        Store posterior predictive samples for each model.

        Parameters
        ----------
        models_dict
            Dictionary of models to store posterior predictive samples for.
        batch_size
            Batch size for generating posterior predictive samples.
        indices
            Indices to generate posterior predictive samples for.
        """
        self.batch_size = batch_size

        samples_dict = {}
        for m, model in self.models.items():
            pp_counts = model.posterior_predictive_sample(
                model.adata,
                n_samples=self.n_samples,
                batch_size=self.batch_size,
                indices=indices,
            )
            if isinstance(pp_counts, dict):
                pp_counts = pp_counts[self.modality]
            samples_dict[m] = DataArray(
                data=pp_counts,
                coords={
                    "cells": list(self.adata.obs_names),
                    "features": list(self.adata.var_names),
                    "samples": np.arange(self.n_samples),
                },
            )
        samples_dict[DATA_VAR_RAW] = DataArray(
            data=self.raw_counts,
            coords={"cells": list(self.adata.obs_names), "features": list(self.adata.var_names)},
        )
        self.samples_dataset = Dataset(samples_dict)

    def coefficient_of_variation(self, dim: Dims = "cells") -> None:
        """
        Calculate the coefficient of variation (CV) for each model and the raw counts.

        The CV is computed over the cells or features dimension per sample. The mean CV is then
        computed over all samples.

        Parameters
        ----------
        dim
            Dimension to compute CV over.
        """
        identifier = METRIC_CV_CELL if dim == "features" else METRIC_CV_GENE
        mean = self.samples_dataset.mean(dim=dim, skipna=False)
        # we use a trick to compute the std to speed it up: std = E[X^2] - E[X]^2
        # a square followed by a sqrt is ok here because this is counts data (no negative values)
        self.samples_dataset = np.square(self.samples_dataset)
        std = np.sqrt(self.samples_dataset.mean(dim=dim, skipna=False) - np.square(mean))
        self.samples_dataset = np.sqrt(self.samples_dataset)
        # now compute the CV
        cv = std / mean
        # It's ok to make things dense here
        cv = _make_dataset_dense(cv)
        cv_mean = cv.mean(dim="samples", skipna=True)
        cv_mean[DATA_VAR_RAW].data = np.nan_to_num(cv_mean[DATA_VAR_RAW].data)
        self.metrics[identifier] = cv_mean.to_dataframe()

    def zero_fraction(self) -> None:
        """Fraction of zeros in raw counts for a specific gene"""
        pp_samples = self.samples_dataset
        mean = (pp_samples != 0).mean(dim="cells", skipna=False).mean(dim="samples", skipna=False)
        mean = _make_dataset_dense(mean)
        self.metrics[METRIC_ZERO_FRACTION] = mean.to_dataframe()

    def calibration_error(self, confidence_intervals: list[float] | float = None) -> None:
        """Calibration error for each observed count.

        For a series of credible intervals of the samples, the fraction of observed counts that
        fall within the credible interval is computed. The calibration error is then the squared
        difference between the observed fraction and the true interval width.

        For this metric, lower is better.

        Parameters
        ----------
        confidence_intervals
            List of confidence intervals to compute calibration error for.
            E.g., [0.01, 0.02, 0.98, 0.99]

        Notes
        -----
        This does not work on sparse data and can cause large memory usage.
        """
        if confidence_intervals is None:
            ps = [2.5, 5, 7.5, 10, 12.5, 15, 17.5, 82.5, 85, 87.5, 90, 92.5, 95, 97.5]
            ps = [p / 100 for p in ps]
        else:
            if len(confidence_intervals) % 2 != 0:
                raise ValueError("Confidence intervals must be even")
            ps = confidence_intervals
        pp_samples = self.samples_dataset
        # TODO: Reimplement to work on sparse data
        pp_samples = _make_dataset_dense(pp_samples)
        # results in (quantiles, cells, features)
        quants = pp_samples.quantile(q=ps, dim="samples", skipna=False)
        credible_interval_indices = [(i, len(ps) - (i + 1)) for i in range(len(ps) // 2)]

        model_cal = {}
        for model in pp_samples.data_vars:
            if model == DATA_VAR_RAW:
                continue
            cal_error_features = 0
            for interval in credible_interval_indices:
                start = interval[0]
                end = interval[1]
                true_width = ps[end] - ps[start]
                greater_than = (quants[DATA_VAR_RAW] >= quants.model1.isel(quantile=start)).data
                less_than = (quants[DATA_VAR_RAW] <= quants.model1.isel(quantile=end)).data
                # Logical and
                ci = greater_than * less_than
                pci_features = ci.mean()
                cal_error_features += (pci_features - true_width) ** 2
            model_cal[model] = {
                "features": cal_error_features,
            }
        self.metrics[METRIC_CALIBRATION] = pd.DataFrame.from_dict(model_cal)

    @dependencies("scanpy")
    def differential_expression(
        self,
        de_groupby: str | None = None,
        de_method: str = "t-test",
        n_samples: int = 1,
        cell_scale_factor: float = 1e4,
        p_val_thresh: float = DEFAULT_DE_P_VAL_THRESHOLD,
        n_top_genes_fallback: int = DEFAULT_DE_N_TOP_GENES_OVERLAP,
    ):
        """
        Compute differential expression (DE) metrics.

        If n_samples > 1, all metrics are averaged over a posterior predictive dataset.

        Parameters
        ----------
        de_groupby
            The column name in `adata_obs_raw` that contains the groupby information.
        de_method
            The DE method to use. See :meth:`~scanpy.tl.rank_genes_groups` for more details.
        n_samples
            The number of posterior predictive samples to use for the DE analysis.
        cell_scale_factor
            The cell scale factor to use for normalization before DE.
        p_val_thresh
            The p-value threshold to use for the DE analysis.
        n_top_genes_fallback
            The number of top genes to use for the DE analysis if the number of genes
            with a p-value < p_val_thresh is zero.
        """
        if 10 * n_top_genes_fallback > self.adata.n_vars:
            warnings.warn(
                f"n_top_genes_fallback={n_top_genes_fallback} is greater than 10% of n_vars"
                f" {self.adata.n_vars} in the dataset. Setting it to 10% of n_vars.",
                UserWarning,
                stacklevel=2,
            )
            n_top_genes_fallback = int(0.1 * self.adata.n_vars)

        import scanpy as sc

        if n_samples > self.n_samples:
            raise ValueError(
                f"n_samples={n_samples} is greater than the number of samples already recorded "
                f"({self.n_samples})"
            )
        # run DE with the raw counts
        adata_de = AnnData(
            X=self.raw_counts.to_scipy_sparse().tocsr().copy(),
            obs=self.adata.obs,
            var=self.adata.var,
        )
        sc.pp.normalize_total(adata_de, target_sum=cell_scale_factor)
        sc.pp.log1p(adata_de)
        if de_groupby is None:
            sc.tl.pca(adata_de)
            sc.pp.neighbors(adata_de)
            sc.tl.leiden(adata_de, key_added="leiden_scvi_criticism")
            de_groupby = "leiden_scvi_criticism"
        with warnings.catch_warnings():
            warnings.simplefilter(action="ignore", category=pd.errors.PerformanceWarning)
            sc.tl.rank_genes_groups(
                adata_de,
                de_groupby,
                use_raw=False,
                method=de_method,
                key_added=UNS_NAME_RGG_RAW,
            )

        # get posterior predictive samples from the model (aka approx counts)
        pp_samples = self.samples_dataset
        # create adata object to run DE on the approx counts
        # X here will be overwritten
        adata_approx = AnnData(X=adata_de.X, obs=adata_de.obs, var=adata_de.var)
        de_keys = {}
        models = [model for model in pp_samples.data_vars if model != DATA_VAR_RAW]
        for model in models:
            if model not in de_keys:
                de_keys[model] = []
            for k in range(n_samples):
                one_sample = pp_samples[model].isel(samples=k)
                # overwrite X with the posterior predictive sample
                # This allows us to save all the DE results in the same adata object
                one_sample_data = (
                    one_sample.data.to_scipy_sparse().tocsr()
                    if isinstance(one_sample.data, SparseArray)
                    else one_sample
                )
                adata_approx.X = one_sample_data.copy()
                sc.pp.normalize_total(adata_approx, target_sum=cell_scale_factor)
                sc.pp.log1p(adata_approx)

                # run DE with the imputed normalized data
                with warnings.catch_warnings():
                    warnings.simplefilter(action="ignore", category=pd.errors.PerformanceWarning)
                    key_added = f"{UNS_NAME_RGG_PPC}_{model}_sample_{k}"
                    de_keys[model].append(key_added)
                    sc.tl.rank_genes_groups(
                        adata_approx,
                        de_groupby,
                        use_raw=False,
                        method=de_method,
                        key_added=key_added,
                    )

        groups = adata_de.obs[de_groupby].astype("category").cat.categories
        cell_counts = adata_de.obs[de_groupby].value_counts()
        df = pd.DataFrame(
            index=np.arange(len(groups) * len(models)),
            columns=[
                "gene_f1",
                "lfc_mae",
                "lfc_pearson",
                "lfc_spearman",
                "roc_auc",
                "pr_auc",
                "group",
                "model",
                "n_cells",
            ],
        )
        # Initialize storage for metrics
        self.metrics[METRIC_DIFF_EXP] = {"lfc_per_model_per_group": {}}

        for i, group in enumerate(groups):
            raw_group_data = sc.get.rank_genes_groups_df(
                adata_de, group=group, key=UNS_NAME_RGG_RAW
            )
            raw_group_data.set_index("names", inplace=True)

            for model, model_keys in de_keys.items():
                # Storage for metrics across samples
                gene_overlap_f1s, lfc_maes, lfc_pearsons, lfc_spearmans = [], [], [], []
                roc_aucs, pr_aucs, rgds, sgds = [], [], [], []
                for de_key in model_keys:
                    sample_group_data = sc.get.rank_genes_groups_df(
                        adata_approx, group=group, key=de_key
                    )
                    sample_group_data.set_index("names", inplace=True)

                    # Gene Overlap F1
                    all_genes = raw_group_data.index
                    top_genes_raw = raw_group_data[:n_top_genes_fallback].index
                    top_genes_sample = sample_group_data[:n_top_genes_fallback].index
                    true_genes = np.isin(all_genes, top_genes_raw).astype(int)
                    pred_genes = np.isin(all_genes, top_genes_sample).astype(int)
                    gene_overlap_f1s.append(
                        precision_recall_fscore_support(true_genes, pred_genes, average="binary")[
                            2
                        ]
                    )
                    # Log-fold change (LFC) metrics
                    sample_group_data = sample_group_data.reindex(raw_group_data.index)
                    rgd, sgd = (
                        raw_group_data["logfoldchanges"],
                        sample_group_data["logfoldchanges"],
                    )
                    rgds.append(rgd)
                    sgds.append(sgd)
                    lfc_maes.append(np.mean(np.abs(rgd - sgd)))
                    lfc_pearsons.append(pearsonr(rgd, sgd)[0])
                    lfc_spearmans.append(spearmanr(rgd, sgd)[0])

                    # ROC and PR metrics
                    raw_adj_p_vals = raw_group_data["pvals_adj"]
                    true = raw_adj_p_vals < p_val_thresh
                    pred = sample_group_data["scores"]

                    # Fallback for no true DE genes and most genes DE.
                    if true.sum() == 0 or true.sum() > (0.5 * len(true)):
                        true = np.zeros_like(pred)
                        true[np.argsort(raw_adj_p_vals)[:n_top_genes_fallback]] = 1

                    roc_aucs.append(roc_auc_score(true, pred))
                    pr_aucs.append(average_precision_score(true, pred))

                # Compute means over samples
                df.loc[i, "model"] = model
                df.loc[i, "group"] = group
                df.loc[i, "gene_f1"] = np.mean(gene_overlap_f1s)
                df.loc[i, "lfc_mae"] = np.mean(lfc_maes)
                df.loc[i, "lfc_pearson"] = np.mean(lfc_pearsons)
                df.loc[i, "lfc_spearman"] = np.mean(lfc_spearmans)
                df.loc[i, "roc_auc"] = np.mean(roc_aucs)
                df.loc[i, "pr_auc"] = np.mean(pr_aucs)
                df.loc[i, "n_cells"] = cell_counts[group]

                # Store LFCs for raw vs approx
                rgd_avg, sgd_avg = pd.DataFrame(rgds).mean(axis=0), pd.DataFrame(sgds).mean(axis=0)
                self.metrics[METRIC_DIFF_EXP]["lfc_per_model_per_group"].setdefault(model, {})
                self.metrics[METRIC_DIFF_EXP]["lfc_per_model_per_group"][model][str(group)] = (
                    pd.DataFrame([rgd_avg, sgd_avg], index=["raw", "approx"]).T
                )

        self.metrics[METRIC_DIFF_EXP]["summary"] = df
