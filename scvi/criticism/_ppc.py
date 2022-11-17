import logging
from typing import Dict, Optional, Sequence, Type, Union

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from anndata import AnnData
from scipy.sparse import coo_matrix, csr_matrix
from scipy.stats import mannwhitneyu, pearsonr, spearmanr
from sklearn.metrics import mean_absolute_error as mae
from sklearn.metrics import mean_squared_error as mse

import scvi
from scvi.model.base import BaseModelClass

logger = logging.getLogger(__name__)


class PPC:
    """Posterior predictive checks for comparing single-cell generative models"""

    def __init__(
        self,
        n_samples: int = 1,
        raw_counts: Optional[Union[np.ndarray, csr_matrix, coo_matrix]] = None,
    ):
        if isinstance(raw_counts, np.ndarray):
            self.raw_counts = coo_matrix(raw_counts)
        elif isinstance(raw_counts, csr_matrix):
            self.raw_counts = raw_counts.tocoo()
        elif isinstance(raw_counts, coo_matrix):
            self.raw_counts = raw_counts
        else:
            self.raw_counts = None
        self.posterior_predictive_samples = {}
        self.n_samples = n_samples
        self.models = {}
        self.metrics = {}

    def store_scvi_posterior_samples(
        self,
        models_dict: Dict[str, Type[BaseModelClass]],
        batch_size=32,
        indices=None,
    ):
        """Gathers posterior predictive samples."""
        self.models = models_dict
        self.batch_size = batch_size
        first_model = next(iter(models_dict.keys()))
        self.dataset = models_dict[first_model].adata

        for m, model in self.models.items():
            pp_counts = model.posterior_predictive_sample(
                model.adata,
                n_samples=self.n_samples,
                batch_size=self.batch_size,
                indices=indices,
            )
            self.posterior_predictive_samples[m] = pp_counts

    def coefficient_of_variation(self, cell_wise: bool = True):
        """
        Calculate the coefficient of variation.

        Parameters:
            cell_wise: Calculate for each cell across genes if True, else do the reverse.
        """
        axis = 1 if cell_wise is True else 0
        identifier = "cv_cell" if cell_wise is True else "cv_gene"
        df = pd.DataFrame()
        pp_samples = self.posterior_predictive_samples.items()
        for m, samples in pp_samples:
            cv = np.nanmean(
                np.std(samples, axis=axis) / np.mean(samples, axis=axis),
                axis=-1,
            )

            df[m] = cv.ravel()
            df[m] = np.nan_to_num(df[m])

        raw = self.raw_counts.todense()
        df["Raw"] = pd.DataFrame(
            np.asarray(np.std(raw, axis=axis)).squeeze()
            / np.asarray(np.mean(raw, axis=axis)).squeeze()
        )
        df["Raw"] = np.nan_to_num(df["Raw"])

        self.metrics[identifier] = df

    def mann_whitney_u(self):
        """Calculate the Mann–Whitney U statistic."""
        feat_df = pd.DataFrame()
        pp_samples = self.posterior_predictive_samples.items()
        raw = self.raw_counts.todense()
        for m, samples in pp_samples:
            sam = samples
            feats = []
            for g in range(samples.shape[1]):
                Us = []
                for n in range(samples.shape[2]):
                    U, _ = mannwhitneyu(sam[:, g, n], raw[:, g])
                    Us.append(U)
                feats.append(np.mean(Us))
            to_add = feats
            if len(to_add) != raw.shape[1]:
                raise ValueError()
            feat_df[m] = to_add
        self.metrics["mannwhitneyu"] = feat_df


def _get_ppc_metrics(
    adata: AnnData,
    model: Type[BaseModelClass],
    metric: str,
    n_samples: int,
    indices: Sequence[int],
    layer: Optional[str] = None,
) -> None:
    raw_data = adata[indices, :].X if layer is None else adata[indices, :].layers[layer]

    sp = PPC(n_samples=n_samples, raw_counts=raw_data)
    model_name = f"{model.__class__.__name__}"
    models_dict = {model_name: model}
    sp.store_scvi_posterior_samples(models_dict, indices=indices)

    model_metric, raw_metric = None, None
    if metric == "cv_cell" or metric == "cv_gene":
        sp.coefficient_of_variation(cell_wise=(metric == "cv_cell"))
        model_metric = sp.metrics[metric][model_name].values
        raw_metric = sp.metrics[metric]["Raw"].values
    elif metric == "mannwhitneyu":
        sp.mann_whitney_u()
        model_metric = sp.metrics[metric][model_name].values
        raw_metric = None
    else:
        raise NotImplementedError(f"Unknown metric: {metric}")

    return sp, model_metric, raw_metric


def plot_ppc(title: str, model_metric, raw_metric, metric: str):
    """Plot and log summary ppc results for the given model and raw metric vectors"""
    # from https://stackoverflow.com/a/28216751
    def add_identity(axes, *line_args, **line_kwargs):
        (identity,) = axes.plot([], [], *line_args, **line_kwargs)

        def callback(axes):
            low_x, high_x = axes.get_xlim()
            low_y, high_y = axes.get_ylim()
            low = max(low_x, low_y)
            high = min(high_x, high_y)
            identity.set_data([low, high], [low, high])

        callback(axes)
        axes.callbacks.connect("xlim_changed", callback)
        axes.callbacks.connect("ylim_changed", callback)
        return axes

    if metric == "cv_cell" or metric == "cv_gene":
        # mae, mse, pearson corr, spearman corr
        logger.info(
            f"{title}:\n"
            f"Mean Absolute Error={mae(model_metric, raw_metric):.2f},\n"
            f"Mean Squared Error={mse(model_metric, raw_metric):.2f}\n"
            f"Pearson correlation={pearsonr(model_metric, raw_metric)[0]:.2f}\n"
            f"Spearman correlation={spearmanr(model_metric, raw_metric)[0]:.2f}\n"
        )
        # visual correlation
        plt.scatter(model_metric, raw_metric)
        ax = plt.gca()
        add_identity(ax, color="r", ls="--", alpha=0.5)
        plt.xlabel("model")
        plt.ylabel("raw")
        plt.title(title)
        plt.show()
    elif metric == "mannwhitneyu":
        _, ax = plt.subplots(2, 1, figsize=(10, 12.5), sharex=False)
        sns.boxplot(
            data=np.log10(model_metric),
        )
    else:
        raise NotImplementedError()


def run_ppc(
    adata: AnnData,
    model: Type[BaseModelClass],
    metric: str,
    n_samples: int,
    custom_indices: Optional[Sequence[int]] = None,
    n_indices: Optional[int] = None,
    layer: Optional[str] = None,
    do_plot: bool = True,
):
    """Compute the given PPC metric for the given model, data and indices. Plot results by default"""
    if scvi.data._utils._get_latent_adata_type(adata) is not None:
        raise ValueError("Please provide the anndata object containing the full counts")

    # determine indices to use
    if custom_indices is not None:
        indices = custom_indices
    else:
        if n_indices is not None:
            indices = np.random.randint(0, adata.n_obs, n_indices)
        else:
            indices = np.arange(adata.n_obs)

    sp, model_metrics, raw_metrics = _get_ppc_metrics(
        adata,
        model,
        metric,
        n_samples=n_samples,
        indices=indices,
        layer=layer,
    )

    if do_plot:
        model_name = f"{model.__class__.__name__}"
        plot_ppc(
            f"model={model_name} | metric={metric} | n_cells={len(indices)}",
            model_metrics,
            raw_metrics,
            metric,
        )

    return sp
