import logging
from typing import Union

import numpy as np
import pandas as pd
import scipy
import torch
from anndata import AnnData
from scipy.optimize import linear_sum_assignment
from sklearn.neighbors import NearestNeighbors

logger = logging.getLogger(__name__)


def nearest_neighbor_overlap(x1, x2, k=100):
    """
    Compute the overlap between the k-nearest neighbor graph of x1 and x2.

    Using Spearman correlation of the adjacency matrices.
    Compute the overlap fold enrichment between the protein and mRNA-based cell 100-nearest neighbor
        graph and the Spearman correlation of the adjacency matrices.
    """
    if len(x1) != len(x2):
        raise ValueError("len(x1) != len(x2)")
    n_samples = len(x1)
    k = min(k, n_samples - 1)
    nne = NearestNeighbors(n_neighbors=k + 1)  # "n_jobs=8
    nne.fit(x1)
    kmatrix_1 = nne.kneighbors_graph(x1) - scipy.sparse.identity(n_samples)
    nne.fit(x2)
    kmatrix_2 = nne.kneighbors_graph(x2) - scipy.sparse.identity(n_samples)

    # 1 - spearman correlation from knn graphs
    spearman_correlation = scipy.stats.spearmanr(
        kmatrix_1.A.flatten(), kmatrix_2.A.flatten()
    )[0]
    # 2 - fold enrichment
    set_1 = set(np.where(kmatrix_1.A.flatten() == 1)[0])
    set_2 = set(np.where(kmatrix_2.A.flatten() == 1)[0])
    fold_enrichment = (
        len(set_1.intersection(set_2))
        * n_samples**2
        / (float(len(set_1)) * len(set_2))
    )
    return spearman_correlation, fold_enrichment


def unsupervised_clustering_accuracy(
    y: Union[np.ndarray, torch.Tensor], y_pred: Union[np.ndarray, torch.Tensor]
) -> tuple:
    """Unsupervised Clustering Accuracy."""
    if len(y_pred) != len(y):
        raise ValueError("len(y_pred) != len(y)")
    u = np.unique(np.concatenate((y, y_pred)))
    n_clusters = len(u)
    mapping = dict(zip(u, range(n_clusters)))
    reward_matrix = np.zeros((n_clusters, n_clusters), dtype=np.int64)
    for y_pred_, y_ in zip(y_pred, y):
        if y_ in mapping:
            reward_matrix[mapping[y_pred_], mapping[y_]] += 1
    cost_matrix = reward_matrix.max() - reward_matrix
    row_assign, col_assign = linear_sum_assignment(cost_matrix)

    # Construct optimal assignments matrix
    row_assign = row_assign.reshape((-1, 1))  # (n,) to (n, 1) reshape
    col_assign = col_assign.reshape((-1, 1))  # (n,) to (n, 1) reshape
    assignments = np.concatenate((row_assign, col_assign), axis=1)

    optimal_reward = reward_matrix[row_assign, col_assign].sum() * 1.0
    return optimal_reward / y_pred.size, assignments


def knn_purity(latent, label, n_neighbors=30):
    nbrs = NearestNeighbors(n_neighbors=n_neighbors + 1).fit(latent)
    indices = nbrs.kneighbors(latent, return_distance=False)[:, 1:]
    neighbors_labels = np.vectorize(lambda i: label[i])(indices)

    # pre cell purity scores
    scores = ((neighbors_labels - label.reshape(-1, 1)) == 0).mean(axis=1)
    res = [
        np.mean(scores[label == i]) for i in np.unique(label)
    ]  # per cell-type purity

    return np.mean(res)

def silhouette_labels(
    adata: AnnData,
    labels_key: str = "labels_bench",
    metric: str = "euclidean",
    embedding: str = "X_scvi",
    scale: bool = True,
) -> float:
    """
    Wrapper of Silhoutte score from sklearn, with respect to observed cell type labels.
    The score is scaled into the range [0, 1], where larger values indicate better, denser clusters.
    """
    from sklearn.metrics import silhouette_score

    if embedding not in adata.obsm.keys():
        raise KeyError(f"{embedding} not in obsm")
    asw = silhouette_score(adata.obsm[embedding], adata.obs[labels_key], metric=metric)
    if scale:
        asw += 1
        asw /= 2
    return asw

def silhouette_batch(
    adata: AnnData,
    batch_key: str = "batch_bench",
    labels_key: str = "labels_bench",
    metric: str = "euclidean",
    embedding: str = "X_scvi",
    scale: bool = True,
):
    """
    Silhouette score of batch labels subsetted for each group.
    batch_ASW = 1 - abs(ASW), so score of 1 means ideally mixed, 0, not mixed.
    Parameters
    ----------
    adata
        AnnData object
    batch_key
        Key in `adata.obs` containing batch information
    labels_key
        Key in `adata.obs` containing cell type labels
    metric
        A valid sklearn metric
    embedding
        Key in `adata.obsm` containing scvi-tools generated low-dim representation
    Returns
    -------
    all_scores
        Absolute silhouette scores per group label
    group_means
        Mean per cell type (group)
    """
    from sklearn.metrics import silhouette_score

    if embedding not in adata.obsm.keys():
        raise KeyError(f"{embedding} not in obsm")

    groups = adata.obs[labels_key].unique()
    sil_all = pd.DataFrame(columns=["silhouette_score"], index=groups)

    for group in groups:
        adata_group = adata[adata.obs[labels_key] == group]
        if adata_group.obs[batch_key].nunique() == 1:
            continue
        sil_per_group = silhouette_score(
            adata_group.obsm[embedding], adata_group.obs[batch_key], metric=metric
        )
        if scale:
            sil_per_group = (sil_per_group + 1) / 2
            sil_per_group = np.abs(sil_per_group)
            sil_per_group = 1 - sil_per_group

        sil_all.loc[group, "silhouette_score"] = sil_per_group
    return sil_all

def lisi(
    adata: AnnData,
    batch_key: str = "batch_bench",
    labels_key: str = "labels_bench",
    embedding: str = "X_scvi",
    n_jobs: int = 1,
) -> float:
    """
    Wrapper of compute_lisi from harmonypy package. Higher is better.
    Suppose one of the columns in metadata is a categorical variable with 3 categories.
        - If LISI is approximately equal to 3 for an item in the data matrix,
          that means that the item is surrounded by neighbors from all 3
          categories.
        - If LISI is approximately equal to 1, then the item is surrounded by
          neighbors from 1 category.
    Returns
    -------
    iLISI
        lisi computed w.r.t. batches
    cLISI
        lisi computed w.r.t. labels
    """
    from harmonypy import compute_lisi
    from sklearn.utils import parallel_backend

    if embedding not in adata.obsm.keys():
        raise KeyError(f"{embedding} not in obsm")

    if n_jobs != 1:
        with parallel_backend("threading", n_jobs=n_jobs):
            lisi = compute_lisi(
                adata.obsm[embedding],
                metadata=adata.obs,
                label_colnames=[batch_key, labels_key],
            )
    else:
        lisi = compute_lisi(
            adata.obsm[embedding],
            metadata=adata.obs,
            label_colnames=[batch_key, labels_key],
        )
    return lisi

def compute_ari(adata: AnnData, labels_key: str = "labels_bench") -> float:
    """
    Adjusted rand index computed at various clustering resolutions, max reported.
    Parameters
    ----------
    adata
        AnnData object
    labels_key
        Key in `adata.obs` containing cell type labels
    """
    from sklearn.metrics.cluster import adjusted_rand_score

    group1 = adata.obs[labels_key]
    aris = []
    keys = adata.obs.keys()
    for k in keys:
        if "leiden_bench" in k:
            group2 = adata.obs[k]
            ari = adjusted_rand_score(group1, group2)
            aris.append(ari)

    return np.max(aris)

def compute_nmi(adata: AnnData, labels_key: str = "labels_bench") -> float:
    """
    Adjusted rand index computed at various clustering resolutions, max reported.
    Parameters
    ----------
    adata
        AnnData object
    labels_key
        Key in `adata.obs` containing cell type labels
    """
    from sklearn.metrics.cluster import normalized_mutual_info_score

    group1 = adata.obs[labels_key]
    nmis = []
    keys = adata.obs.keys()
    for k in keys:
        if "leiden_bench" in k:
            group2 = adata.obs[k]
            nmi = normalized_mutual_info_score(group1, group2)
            nmis.append(nmi)

    return np.max(nmis)
