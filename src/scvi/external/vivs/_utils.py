from __future__ import annotations

from typing import TYPE_CHECKING

import numpy as np
import scanpy as sc
from sklearn.cluster import KMeans
from sklearn.decomposition import PCA
from sklearn.model_selection import ParameterGrid

if TYPE_CHECKING:
    from anndata import AnnData


def select_genes(
    adata: AnnData,
    n_top_genes: int,
    preselected_genes: list[str] | None = None,
    seed: int = 0,
) -> AnnData:
    """Select a representative gene subset via a highly-variable-genes + clustering heuristic.

    Ported from VIVS's original JAX implementation (:cite:p:`BoyeauVIVS24`); pure
    scanpy/sklearn, unchanged by the torch port.

    Parameters
    ----------
    adata
        AnnData with raw counts in ``.X``.
    n_top_genes
        Number of genes to keep.
    preselected_genes
        Gene names to always keep in addition to the selected ones.
    seed
        Random seed for the KMeans clustering step.
    """
    adata_ = adata.copy()
    preselected_genes = preselected_genes if preselected_genes is not None else []

    adata_log = adata_.copy()
    sc.pp.normalize_total(adata_log, target_sum=1e6)
    sc.pp.log1p(adata_log)
    pca_ = PCA(n_components=50).fit(adata_log.X)
    sc.pp.highly_variable_genes(adata_, n_top_genes=n_top_genes, flavor="seurat_v3")

    clusters = KMeans(n_clusters=n_top_genes, random_state=seed, n_init=10).fit_predict(
        pca_.components_.T
    )
    adata_.var.loc[:, "clusters"] = clusters
    adata_.var.index.name = "index"
    selected_genes = (
        adata_.var.reset_index()
        .groupby("clusters")
        .apply(lambda x: x.sort_values("variances_norm").iloc[-1]["index"])
        .values
    )
    union_genes = np.union1d(selected_genes, preselected_genes)
    return adata_[:, adata_.var.index.isin(union_genes)].copy()


def select_architecture(
    adata: AnnData,
    y_obsm_key: str,
    xy_model_kwargs_grid: dict,
    batch_key: str | None = None,
    max_epochs: int = 1,
    **vivs_kwargs,
) -> dict:
    """Grid search over the importance-score net's architecture, by validation loss.

    Trains only the ``xy`` phase for each candidate (each candidate gets its own fresh
    generative VAE, since VIVS's constructor always builds one unless an ``x_model`` is
    passed) and picks the combination with the lowest final training-loss.
    """
    from scvi.external.vivs._model import VIVS

    parameter_grid = list(ParameterGrid(xy_model_kwargs_grid))
    results = []
    for grid_params in parameter_grid:
        VIVS.setup_anndata(adata, y_obsm_key=y_obsm_key, batch_key=batch_key)
        model = VIVS(adata, xy_model_kwargs=grid_params, **vivs_kwargs)
        # check_val_every_n_epoch=1 is required: scvi-tools' Trainer does not validate every
        # epoch by default (confirmed empirically — with max_epochs=1 and no override,
        # model.history_ contains no "elbo_validation" key at all), so a small max_epochs
        # grid search would silently have no validation loss to compare without this.
        model.train(max_epochs=max_epochs, check_val_every_n_epoch=1)
        # "elbo_validation" is the key TrainingPlan.compute_and_log_metrics always logs
        # (mode="validation"), regardless of what a module's loss() actually optimizes —
        # confirmed via scvi.train.TrainingPlan.compute_and_log_metrics source.
        final_val_loss = model.history_["elbo_validation"].iloc[-1, 0]
        results.append({**grid_params, "val_loss": final_val_loss})

    best = min(results, key=lambda r: r["val_loss"])
    return {k: v for k, v in best.items() if k != "val_loss"}
